# Hi-C registered protocol

A Nextflow workflow (`main.nf`) that turns raw alignment-pair files into `.pairs` files and per-library QC, and a delivery
script (`deliver.sh`) that runs it and reports the analysis parameters.

## Run

```
bash deliver.sh --samplesheet SHEET.csv --annotations DIR --out OUT [--settings JSON] [--assay insitu|capture]
                [--baitmap FILE] [--keep-rejected]
```

Environment: `NEXTFLOW` is the workflow engine launcher (`nextflow` if unset), `PYTHON` the interpreter (`python3` if
unset); the analysis package `hicparams` must be importable by it.

| Option | Meaning |
|---|---|
| `--samplesheet` | CSV with columns `library`, `condition`, `raw`; `raw` is a path to the raw file, relative to the samplesheet |
| `--annotations` | directory with `fragments.bed`, `chrom.sizes`, `ep_pairs.bedpe`, `domains.bed` |
| `--out` | output directory |
| `--settings` | assay settings (`assay_settings.json` next to `deliver.sh` if omitted) |
| `--assay` | `insitu` (if omitted) or `capture`; `capture` needs `--baitmap` |
| `--baitmap` | BED of baited fragments (chromosome, start, end, fragment name, bait name) |
| `--keep-rejected` | also keep the rejected pairs of each class for every library |

Examples:

```
bash deliver.sh --samplesheet inputs/samplesheet.csv --annotations inputs/annotations --out /tmp/run1
bash deliver.sh --samplesheet inputs/samplesheet.csv --annotations inputs/annotations --out /tmp/run2 \
    --assay capture --baitmap inputs/designs/panel.baitmap
bash deliver.sh --samplesheet inputs/samplesheet.csv --annotations inputs/annotations --out /tmp/run3 --keep-rejected
```

`main.nf` can also be run directly with `nextflow run main.nf --samplesheet ... --settings ... --annotations ...
[--assay ...] [--baitmap ...] [--keep_rejected true]`.

## Inputs

* Raw file: tab separated, no header, columns `read_id chrom1 pos1 strand1 mapq1 chrom2 pos2 strand2 mapq2`, gzip.
* `assay_settings.json`: settings read by the stage scripts in `bin/`.

## Stages (`bin/`)

`hic_parse.py` (raw file to `.pairs`), `hic_dedup.py`, `hic_filter.py`, `hic_qc.py` (per-library counts and the
library-complexity verdict), `qc_table.py` (QC files to one table). Each prints its usage with `--help`.

## Outputs

Under `OUT/results`:

* `pairs.csv`: one row per library of the samplesheet, columns `library, condition, passed, parsed, dedup, filtered,
  classes, design`.
* `qc.csv`: one row per library, columns `library, condition, passed, qc`.
* `pairs/parsed/`, `pairs/dedup/`, `pairs/filtered/`, `classes/`, `qc/`: the files the rows name.

Under `OUT`: `parameters.tsv` (one row per library of `pairs.csv`, from `python -m hicparams`), `qc.tsv`, `run.log`.
