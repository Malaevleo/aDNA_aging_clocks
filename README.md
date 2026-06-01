# Aging Clocks prediction for aDNA BED methylation maps
<img width="1920" height="1080" alt="gitpic" src="https://github.com/user-attachments/assets/e9ff233f-641c-4d0a-9663-ce2d7973c286" />

This repository contains a standalone script, `predict_age_from_bed.py`, that:

1. converts one or more BED files with columns `chr`, `pos`, and `beta` into Illumina-style CpG tables. Pay special attention to the format of BED files.
2. predicts biological age with `pyaging`
3. writes per-sample CSV outputs and, for multiple inputs, one combined results table

## Requirements

Install the pinned dependencies from `requirements.txt`:

```powershell
pip install -r requirements.txt
```

## EPIC Manifest Used For Conversion

The BED to Illumina conversion uses this EPIC manifest path by default:

```text
./GSE138307/IDAT/GPL21145_MethylationEPIC_15073387_v-1-0.csv
```

You can override it with `--manifest` if needed. We suggest using the one that comes from GSE138307. It can be found here: [https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138307](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138307)

## Single-File Usage

```powershell
python predict_age_from_bed.py ./sample.dammet.bed
```

This writes:

- `sample_illumina.csv`
- `sample_pyaging_predictions.csv`

## Reuse Existing Illumina CSV

```powershell
python predict_age_from_bed.py ./sample.dammet.bed --reuse-illumina --illumina-out ./sample_illumina.csv
```

## Multiple BED Files

```powershell
python predict_age_from_bed.py ./sample.dammet.bed ./sample.roam.bed --combined-out ./combined_pyaging_predictions.csv
```

This writes the usual per-sample outputs and also a combined CSV with all predictions.

## Output Columns

Each predictions CSV contains:

- `sample_id`
- `bed_file`
- `epic_manifest`
- `illumina_file`
- `clock`
- `predicted_age`

## Note

- The script runs three clocks: `Horvath2013`, `AltumAge`, and `ZhangEn`. They proved to be most accurate for bone tissue samples
- The script can be applied to any BED file but was intially developed for aDNA studies 
