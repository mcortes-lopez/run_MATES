#!/usr/bin/env python3
"""
Standalone long-read MATES pipeline with optional corrected counting.

Workflow:
1) split BAMs for long-read processing
2) count TE locus reads
   - default MATES counter, or
   - correction logic from an external script
3) normalize corrected column names for MATES compatibility
4) build long-read TE matrix output
"""

from __future__ import annotations

import argparse
import importlib
import importlib.util
import multiprocessing as mp
import os
import sys
import types
from pathlib import Path
from typing import Iterable, List, Tuple


def read_lines(path: Path) -> List[str]:
    with path.open("r", encoding="utf-8") as handle:
        return [line.strip() for line in handle if line.strip()]


def resolve_input_path(raw_path: str, launch_dir: Path) -> Path:
    path = Path(raw_path)
    if not path.is_absolute():
        path = (launch_dir / path).resolve()
    return path


def resolve_ref_paths(te_mode: str, ref_path: str, workdir: Path) -> Tuple[Path, Path]:
    if ref_path == "Default":
        if te_mode == "exclusive":
            csv_path = workdir / "TE_nooverlap.csv"
            bed_path = workdir / "TE_nooverlap.bed"
        else:
            csv_path = workdir / "TE_full.csv"
            bed_path = workdir / "TE_full.bed"
    else:
        user_ref = Path(ref_path)
        if not user_ref.is_absolute():
            user_ref = (workdir / user_ref).resolve()
        suffix = user_ref.suffix.lower()
        if suffix == ".csv":
            csv_path = user_ref
            bed_path = user_ref.with_suffix(".bed")
        elif suffix == ".bed":
            bed_path = user_ref
            csv_path = user_ref.with_suffix(".csv")
        else:
            raise ValueError("`--ref-path` must point to a .csv or .bed file.")

    if not csv_path.exists():
        raise FileNotFoundError(f"Missing TE csv reference: {csv_path}")
    if not bed_path.exists():
        raise FileNotFoundError(f"Missing TE bed reference: {bed_path}")
    return csv_path, bed_path


def load_correction_module(script_path: Path):
    spec = importlib.util.spec_from_file_location("mates_count_correction", script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load correction script: {script_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_mates_modules(mates_root: Path):
    """
    Load only the MATES modules required for long-read quantification, without
    executing MATES/__init__.py (which imports MATES_main and torch).
    """
    mates_pkg_dir = mates_root / "MATES"
    if not mates_pkg_dir.exists():
        raise FileNotFoundError(f"Missing MATES package directory: {mates_pkg_dir}")

    # Build a lightweight package object so submodule imports resolve from disk
    # without triggering MATES/__init__.py.
    if "MATES" not in sys.modules:
        pkg = types.ModuleType("MATES")
        pkg.__path__ = [str(mates_pkg_dir)]  # type: ignore[attr-defined]
        pkg.__file__ = str(mates_pkg_dir / "__init__.py")
        sys.modules["MATES"] = pkg

    bam_processor = importlib.import_module("MATES.bam_processor")
    te_quantifier_longread = importlib.import_module("MATES.TE_quantifier_LongRead")
    return te_quantifier_longread, bam_processor


_CORRECTION_MODULE = None
_TE_DF = None


def _init_correction_worker(correction_script_path: str, bed_ref_path: str) -> None:
    global _CORRECTION_MODULE, _TE_DF
    script = Path(correction_script_path)
    bed = Path(bed_ref_path)
    correction_module = load_correction_module(script)
    raw_te = correction_module.pybedtools.BedTool(str(bed))
    te_df = correction_module.create_TE_vec_df(raw_te)
    _CORRECTION_MODULE = correction_module
    _TE_DF = te_df


def count_output_complete(output_dir: Path) -> bool:
    required = ("read_weights.csv", "TE_unique_Info.csv", "AS_reads.csv")
    return all((output_dir / name).exists() for name in required)


def _count_task_10x(task: Tuple[str, str, str, bool]) -> str:
    sample, barcode, bam_path, skip_existing = task
    bam = Path(bam_path)
    outdir = Path.cwd() / "count_long_reads" / sample / barcode

    if skip_existing and count_output_complete(outdir):
        return "skipped"
    if not bam.exists():
        return "missing"

    _CORRECTION_MODULE.generate_unique_matric_fast(  # type: ignore[union-attr]
        sample,
        barcode,
        str(bam),
        _TE_DF,  # type: ignore[arg-type]
    )
    return "processed"


def _count_task_smartseq(task: Tuple[str, str, bool]) -> str:
    sample, bam_path, skip_existing = task
    bam = Path(bam_path)
    target = Path.cwd() / "count_long_reads" / sample

    if skip_existing and count_output_complete(target):
        return "skipped"
    if not bam.exists():
        return "missing"

    _CORRECTION_MODULE.generate_unique_matric_fast(  # type: ignore[union-attr]
        sample,
        sample,
        str(bam),
        _TE_DF,  # type: ignore[arg-type]
    )

    nested = Path.cwd() / "count_long_reads" / sample / sample
    target.mkdir(parents=True, exist_ok=True)
    for name in ("TE_unique_Info.csv", "AS_reads.csv", "read_weights.csv"):
        src = nested / name
        if src.exists():
            dst = target / name
            dst.write_bytes(src.read_bytes())
    return "processed"


def _summarize_statuses(statuses: Iterable[str]) -> Tuple[int, int, int]:
    processed = 0
    skipped = 0
    missing = 0
    for status in statuses:
        if status == "processed":
            processed += 1
        elif status == "skipped":
            skipped += 1
        elif status == "missing":
            missing += 1
    return processed, skipped, missing


def iter_te_info_files(data_mode: str, sample_list_file: Path, workdir: Path) -> Iterable[Path]:
    root = workdir / "count_long_reads"
    if not root.exists():
        return []

    if data_mode == "Smart_seq":
        for cell_dir in root.iterdir():
            if cell_dir.is_dir():
                candidate = cell_dir / "TE_unique_Info.csv"
                if candidate.exists():
                    yield candidate
        return

    sample_names = read_lines(sample_list_file)
    for sample in sample_names:
        sample_dir = root / sample
        if not sample_dir.exists():
            continue
        for bc_dir in sample_dir.iterdir():
            if not bc_dir.is_dir():
                continue
            candidate = bc_dir / "TE_unique_Info.csv"
            if candidate.exists():
                yield candidate


def normalize_count_tables(data_mode: str, sample_list_file: Path, workdir: Path) -> int:
    import pandas as pd

    normalized = 0
    for te_file in iter_te_info_files(data_mode, sample_list_file, workdir):
        df = pd.read_csv(te_file)
        renames = {}
        if "te_index" in df.columns and "TE_index" not in df.columns:
            renames["te_index"] = "TE_index"
        if "TE_sense_read_num" in df.columns and "TE_region_read_num" not in df.columns:
            renames["TE_sense_read_num"] = "TE_region_read_num"
        changed = False
        if renames:
            df = df.rename(columns=renames)
            changed = True
        if "TE_index" not in df.columns or "TE_region_read_num" not in df.columns:
            raise ValueError(
                f"{te_file} does not contain required columns after normalization. "
                "Expected TE_index and TE_region_read_num."
            )
        target_cols = ["TE_index", "TE_region_read_num"]
        if list(df.columns) != target_cols:
            changed = True
            df = df[target_cols]
        if changed:
            df[target_cols].to_csv(te_file, index=False)
            normalized += 1
    return normalized


def run_corrected_count_10x(
    sample_list_file: Path,
    bc_path_file: Path,
    correction_script: Path,
    bed_ref_path: Path,
    workdir: Path,
    workers: int,
    chunksize: int,
    skip_existing: bool,
) -> Tuple[int, int, int]:
    sample_names = read_lines(sample_list_file)
    bc_files = read_lines(bc_path_file)
    if len(sample_names) != len(bc_files):
        raise ValueError(
            "sample list and barcode-path list must have the same number of lines."
        )

    tasks: List[Tuple[str, str, str, bool]] = []
    for sample, bc_file in zip(sample_names, bc_files):
        bc_path = Path(bc_file)
        if not bc_path.is_absolute():
            bc_path = (workdir / bc_path).resolve()
        barcodes = read_lines(bc_path)
        for barcode in barcodes:
            bam = workdir / "long_read" / sample / "by_barcode" / f"{barcode}.bam"
            tasks.append((sample, barcode, str(bam), skip_existing))

    if not tasks:
        return 0, 0, 0

    if workers <= 1:
        _init_correction_worker(str(correction_script), str(bed_ref_path))
        statuses = (_count_task_10x(task) for task in tasks)
        return _summarize_statuses(statuses)

    with mp.Pool(
        processes=workers,
        initializer=_init_correction_worker,
        initargs=(str(correction_script), str(bed_ref_path)),
    ) as pool:
        statuses = pool.imap_unordered(_count_task_10x, tasks, chunksize=chunksize)
        return _summarize_statuses(statuses)


def run_corrected_count_smartseq(
    sample_list_file: Path,
    bam_path_file: Path,
    correction_script: Path,
    bed_ref_path: Path,
    workdir: Path,
    workers: int,
    chunksize: int,
    skip_existing: bool,
) -> Tuple[int, int, int]:
    samples = read_lines(sample_list_file)
    bam_paths = read_lines(bam_path_file)
    if len(samples) != len(bam_paths):
        raise ValueError("sample list and bam-path list must have the same number of lines.")

    tasks: List[Tuple[str, str, bool]] = []
    for sample, bam_file in zip(samples, bam_paths):
        bam = Path(bam_file)
        if not bam.is_absolute():
            bam = (workdir / bam).resolve()
        tasks.append((sample, str(bam), skip_existing))

    if not tasks:
        return 0, 0, 0

    if workers <= 1:
        _init_correction_worker(str(correction_script), str(bed_ref_path))
        statuses = (_count_task_smartseq(task) for task in tasks)
        return _summarize_statuses(statuses)

    with mp.Pool(
        processes=workers,
        initializer=_init_correction_worker,
        initargs=(str(correction_script), str(bed_ref_path)),
    ) as pool:
        statuses = pool.imap_unordered(_count_task_smartseq, tasks, chunksize=chunksize)
        return _summarize_statuses(statuses)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run MATES long-read quantification with optional corrected counting."
    )
    parser.add_argument("--te-mode", choices=["inclusive", "exclusive"], required=True)
    parser.add_argument("--data-mode", choices=["10X", "Smart_seq"], required=True)
    parser.add_argument("--sample-list-file", required=True)
    parser.add_argument("--bam-path-file", required=True)
    parser.add_argument("--bc-path-file", default=None)
    parser.add_argument("--bc-ind", default="CB")
    parser.add_argument("--threads-num", type=int, default=6)
    parser.add_argument(
        "--process-num",
        type=int,
        default=4,
        help="Used by MATES BAM splitting (per-sample split workers).",
    )
    parser.add_argument("--ref-path", default="Default")
    parser.add_argument(
        "--mates-root",
        default="/gpfs/commons/groups/landau_lab/mariela/tools/MATES",
        help="Path to MATES repository root containing the MATES package.",
    )
    parser.add_argument(
        "--correction-script",
        default="/gpfs/commons/home/sreddy/MATES/count_ds_reads_updated3.py",
        help="Path to correction script used for long-read counting.",
    )
    parser.add_argument(
        "--workdir",
        default=".",
        help="Run directory where long_read/count_long_reads/long_read_output are created.",
    )
    parser.add_argument("--skip-split", action="store_true")
    parser.add_argument("--skip-count", action="store_true")
    parser.add_argument("--skip-quantify", action="store_true")
    parser.add_argument("--disable-correction", action="store_true")
    parser.add_argument(
        "--start-from-split-bams",
        action="store_true",
        help="Start from already split barcode BAMs in long_read/<sample>/by_barcode.",
    )
    parser.add_argument(
        "--count-workers",
        type=int,
        default=None,
        help="Workers for corrected counting. Default uses --threads-num.",
    )
    parser.add_argument(
        "--count-chunksize",
        type=int,
        default=20,
        help="Task chunk size for multiprocessing during corrected counting.",
    )
    parser.add_argument(
        "--recount-existing",
        action="store_true",
        help="Recompute outputs even when per-cell corrected files already exist.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    launch_dir = Path.cwd()

    workdir = resolve_input_path(args.workdir, launch_dir)
    sample_list_file = resolve_input_path(args.sample_list_file, launch_dir)
    bam_path_file = resolve_input_path(args.bam_path_file, launch_dir)
    bc_path_file = (
        resolve_input_path(args.bc_path_file, launch_dir) if args.bc_path_file else None
    )
    mates_root = resolve_input_path(args.mates_root, launch_dir)
    correction_script = resolve_input_path(args.correction_script, launch_dir)

    if args.data_mode == "10X" and bc_path_file is None:
        raise ValueError("--bc-path-file is required for data-mode 10X.")
    if args.start_from_split_bams:
        args.skip_split = True

    if not sample_list_file.exists():
        raise FileNotFoundError(f"Missing sample list file: {sample_list_file}")
    if not bam_path_file.exists():
        raise FileNotFoundError(f"Missing bam path file: {bam_path_file}")
    if bc_path_file and not bc_path_file.exists():
        raise FileNotFoundError(f"Missing barcode path file: {bc_path_file}")
    if not mates_root.exists():
        raise FileNotFoundError(f"Missing MATES root: {mates_root}")

    TE_quantifier_LongRead, bam_processor = load_mates_modules(mates_root)

    workdir.mkdir(parents=True, exist_ok=True)
    os.chdir(workdir)
    csv_ref_path, bed_ref_path = resolve_ref_paths(args.te_mode, args.ref_path, workdir)

    print(f"Running in workdir: {workdir}")
    print(f"TE references: csv={csv_ref_path} bed={bed_ref_path}")

    if args.skip_split and not args.skip_count and args.data_mode == "10X":
        missing_split_dirs = []
        for sample in read_lines(sample_list_file):
            split_dir = workdir / "long_read" / sample / "by_barcode"
            if not split_dir.exists():
                missing_split_dirs.append(str(split_dir))
        if missing_split_dirs:
            missing_preview = "\n".join(missing_split_dirs[:5])
            raise FileNotFoundError(
                "Requested to skip split step, but split BAM directories are missing. "
                f"First missing paths:\n{missing_preview}"
            )

    if not args.skip_split:
        print("Step 1/3: split BAM files for long-read processing")
        bam_processor.split_bam_files(
            data_mode=args.data_mode,
            threads_num=args.threads_num,
            sample_list_file=str(sample_list_file),
            bam_path_file=str(bam_path_file),
            process_num=args.process_num,
            bc_ind=args.bc_ind,
            long_read=(args.data_mode == "10X"),
            bc_path_file=str(bc_path_file) if bc_path_file else None,
        )

    if not args.skip_count:
        if args.disable_correction:
            print("Step 2/3: counting long reads with native MATES counter")
            bam_processor.count_long_reads(
                TE_mode=args.te_mode,
                data_mode=args.data_mode,
                threads_num=args.threads_num,
                sample_list_file=str(sample_list_file),
                ref_path=str(csv_ref_path),
                bc_path_file=str(bc_path_file) if bc_path_file else None,
            )
        else:
            print("Step 2/3: counting long reads with correction script")
            if not correction_script.exists():
                raise FileNotFoundError(f"Missing correction script: {correction_script}")
            workers = args.count_workers if args.count_workers is not None else args.threads_num
            workers = max(1, workers)
            chunksize = max(1, args.count_chunksize)
            skip_existing = not args.recount_existing
            if args.data_mode == "10X":
                processed, skipped, missing = run_corrected_count_10x(
                    sample_list_file=sample_list_file,
                    bc_path_file=bc_path_file,  # type: ignore[arg-type]
                    correction_script=correction_script,
                    bed_ref_path=bed_ref_path,
                    workdir=workdir,
                    workers=workers,
                    chunksize=chunksize,
                    skip_existing=skip_existing,
                )
                print(
                    f"Corrected counting complete. Processed={processed}, "
                    f"already_done_skipped={skipped}, missing_bam_skipped={missing}"
                )
            else:
                processed, skipped, missing = run_corrected_count_smartseq(
                    sample_list_file=sample_list_file,
                    bam_path_file=bam_path_file,
                    correction_script=correction_script,
                    bed_ref_path=bed_ref_path,
                    workdir=workdir,
                    workers=workers,
                    chunksize=chunksize,
                    skip_existing=skip_existing,
                )
                print(
                    f"Corrected counting complete. Processed={processed}, "
                    f"already_done_skipped={skipped}, missing_bam_skipped={missing}"
                )

        normalized = normalize_count_tables(args.data_mode, sample_list_file, workdir)
        print(f"Normalized TE_unique_Info.csv files rewritten: {normalized}")

    if not args.skip_quantify:
        print("Step 3/3: building long-read locus TE matrix")
        TE_quantifier_LongRead.quantify_locus_TE_MTX(
            TE_mode=args.te_mode,
            data_mode=args.data_mode,
            sample_list_file=str(sample_list_file),
        )
        print("Finished. Output directory: long_read_output/")


if __name__ == "__main__":
    main()
