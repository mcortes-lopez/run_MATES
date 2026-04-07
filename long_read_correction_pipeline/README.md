# Long-Read MATES Pipeline With Correction

This is a standalone wrapper that runs long-read MATES quantification and can replace the default long-read counting step with:

- `/gpfs/commons/home/sreddy/MATES/count_ds_reads_updated3.py`

The wrapper keeps MATES output compatibility by normalizing corrected count columns to:

- `TE_index`
- `TE_region_read_num`

## Script

- `run_mates_longread_correction.py`

## What It Runs

1. `bam_processor.split_bam_files(..., long_read=True)` for `10X` (regular split for `Smart_seq`)
2. corrected counting (or native `bam_processor.count_long_reads` if disabled)
3. `TE_quantifier_LongRead.quantify_locus_TE_MTX(...)`

## Example (10X, corrected counting)

```bash
python general_scripts/run_MATES/long_read_correction_pipeline/run_mates_longread_correction.py \
  --te-mode inclusive \
  --data-mode 10X \
  --sample-list-file general_scripts/run_MATES/run_files/scRNA_samples.txt \
  --bam-path-file general_scripts/run_MATES/run_files/scRNA_bams.txt \
  --bc-path-file general_scripts/run_MATES/run_files/scRNA_barcodes.txt \
  --bc-ind CB \
  --threads-num 12 \
  --count-workers 12 \
  --process-num 4 \
  --workdir general_scripts/run_MATES
```

## Example (Start from split BAMs)

Use this when `long_read/<sample>/by_barcode/*.bam` already exists and you want to skip re-splitting:

```bash
python general_scripts/run_MATES/long_read_correction_pipeline/run_mates_longread_correction.py \
  --te-mode inclusive \
  --data-mode 10X \
  --sample-list-file general_scripts/run_MATES/run_files/scRNA_samples.txt \
  --bam-path-file general_scripts/run_MATES/run_files/scRNA_bams.txt \
  --bc-path-file general_scripts/run_MATES/run_files/scRNA_barcodes.txt \
  --start-from-split-bams \
  --workdir general_scripts/run_MATES
```

## Example (10X, native MATES counting)

```bash
python general_scripts/run_MATES/long_read_correction_pipeline/run_mates_longread_correction.py \
  --te-mode inclusive \
  --data-mode 10X \
  --sample-list-file general_scripts/run_MATES/run_files/scRNA_samples.txt \
  --bam-path-file general_scripts/run_MATES/run_files/scRNA_bams.txt \
  --bc-path-file general_scripts/run_MATES/run_files/scRNA_barcodes.txt \
  --disable-correction \
  --workdir general_scripts/run_MATES
```

## Notes

- Default `--mates-root` is set to `/gpfs/commons/groups/landau_lab/mariela/tools/MATES`.
- Default correction script path is set to `/gpfs/commons/home/sreddy/MATES/count_ds_reads_updated3.py`.
- For `10X`, `--bc-path-file` is required.
- `--workdir` should contain TE references (`TE_full.csv/.bed` or `TE_nooverlap.csv/.bed`) when using `--ref-path Default`.
- Corrected counting is now parallelized via `--count-workers` and `--count-chunksize`.
- Existing per-cell outputs are skipped by default (faster reruns). Use `--recount-existing` to force recomputation.
- This wrapper avoids importing `MATES_main`, so PyTorch is not required for this long-read counting + matrix workflow.
- `pkg_resources` deprecation messages from `pyranges`/`bam_processor` are warnings, not fatal errors.
