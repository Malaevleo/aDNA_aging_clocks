from __future__ import annotations

"""
BED to pyaging workflow.

This script converts one or more BED files with columns `chr`, `pos`, and `beta`
into Illumina-style CpG tables and predicts biological age with pyaging.

Examples
--------
Single BED file:
    python predict_age_from_bed.py ./sample.dammet.bed

Multiple BED files with one combined output table:
    python predict_age_from_bed.py ./sample.dammet.bed ./sample.roam.bed

Reuse an existing Illumina CSV:
    python predict_age_from_bed.py ./sample.dammet.bed --reuse-illumina --illumina-out ./sample_illumina.csv

EPIC manifest for BED -> Illumina conversion:
    ./GSE138307/IDAT/GPL21145_MethylationEPIC_15073387_v-1-0.csv
"""

import argparse
import gc
from pathlib import Path
from typing import Dict, List

import anndata
import numpy as np
import pandas as pd
import pyaging as pya


BASE_DIR = Path("d:/adna")
EPIC_MANIFEST = BASE_DIR / "GSE138307" / "IDAT" / "GPL21145_MethylationEPIC_15073387_v-1-0.csv"
CHROMS = list(range(1, 23))


def load_epic_manifest(manifest_file: Path) -> pd.DataFrame:
    if not manifest_file.exists():
        raise FileNotFoundError(f"EPIC manifest not found: {manifest_file}")

    ann = pd.read_csv(
        manifest_file,
        skiprows=7,
        usecols=["Name", "CHR", "MAPINFO"],
        low_memory=False,
    )

    ann = ann.rename(columns={"Name": "CpG", "CHR": "chr", "MAPINFO": "pos"})
    ann = ann.dropna(subset=["CpG", "chr", "pos"]).copy()
    ann["chr"] = ann["chr"].astype(str).str.replace("chr", "", regex=False)
    ann = ann[ann["chr"].str.fullmatch(r"\d+")]
    ann = ann[ann["chr"].astype(int).isin(CHROMS)]
    ann["pos"] = pd.to_numeric(ann["pos"], errors="coerce")
    ann = ann.dropna(subset=["pos"]).copy()
    ann["pos"] = ann["pos"].astype(np.int64)
    ann["coord_key"] = ann["chr"] + ":" + ann["pos"].astype(str)
    ann = ann.drop_duplicates(subset=["CpG"], keep="first")
    return ann[["CpG", "coord_key"]]


def convert_bed_to_illumina_with_manifest(
    bed_file: Path,
    out_csv: Path,
    sample_name: str,
    coord_to_cpg: Dict[str, str],
    chunk_size: int = 2_000_000,
) -> None:
    cpg_to_beta: Dict[str, float] = {}

    for chunk in pd.read_csv(bed_file, usecols=["chr", "pos", "beta"], chunksize=chunk_size):
        chunk = chunk.dropna(subset=["chr", "pos", "beta"]).copy()
        chunk["chr"] = chunk["chr"].astype(str).str.replace("chr", "", regex=False)
        chunk["chr"] = pd.to_numeric(chunk["chr"], errors="coerce")
        chunk["pos"] = pd.to_numeric(chunk["pos"], errors="coerce")
        chunk = chunk.dropna(subset=["chr", "pos"]).copy()

        if chunk.empty:
            continue

        chunk["chr"] = chunk["chr"].astype(np.int32)
        chunk["pos"] = chunk["pos"].astype(np.int64)
        coord = chunk["chr"].astype(str) + ":" + chunk["pos"].astype(str)
        cpg = coord.map(coord_to_cpg)
        hit_mask = cpg.notna().values

        if hit_mask.any():
            hit_cpg = cpg[hit_mask].values
            hit_beta = pd.to_numeric(chunk.loc[hit_mask, "beta"], errors="coerce").values
            for cg_id, beta_val in zip(hit_cpg, hit_beta):
                if cg_id not in cpg_to_beta:
                    cpg_to_beta[cg_id] = float(beta_val)

        del chunk, coord, cpg, hit_mask
        gc.collect()

    out_df = pd.DataFrame({"CpG": list(cpg_to_beta.keys()), sample_name: list(cpg_to_beta.values())})
    out_df.to_csv(out_csv, index=False)


def predict_one_sample(csv_file: Path, sample_id: str) -> Dict[str, float]:
    df = pd.read_csv(csv_file)
    if "CpG" not in df.columns or df.shape[1] < 2:
        raise ValueError(f"Invalid illumina CSV format: {csv_file}")

    cpg_col = df["CpG"].astype(str)
    beta_vals = pd.to_numeric(df.iloc[:, 1], errors="coerce").astype(np.float32)
    sample_matrix = pd.DataFrame([beta_vals.values], index=[sample_id], columns=cpg_col.values)

    adata = anndata.AnnData(X=sample_matrix.values)
    adata.obs_names = sample_matrix.index
    adata.var_names = sample_matrix.columns

    aliases = {
        "Horvath2013": ["Horvath2013", "horvath2013"],
        "AltumAge": ["AltumAge", "altumage"],
        "ZhangEn": ["ZhangEn", "zhangen", "ZhangEN"],
    }
    results: Dict[str, float] = {}

    for canonical_clock, alias_list in aliases.items():
        pred_value = np.nan
        last_error = None

        for alias in alias_list:
            try:
                tmp = adata.copy()
                pya.pred.predict_age(tmp, [alias])
                cols = list(tmp.obs.columns)
                cols_lu = {column.lower(): column for column in cols}

                if alias.lower() in cols_lu:
                    pred_value = float(tmp.obs[cols_lu[alias.lower()]].iloc[0])
                    del tmp
                    break
                if canonical_clock.lower() in cols_lu:
                    pred_value = float(tmp.obs[cols_lu[canonical_clock.lower()]].iloc[0])
                    del tmp
                    break

                numeric_cols = [column for column in cols if pd.api.types.is_numeric_dtype(tmp.obs[column])]
                if numeric_cols:
                    pred_value = float(tmp.obs[numeric_cols[0]].iloc[0])
                    del tmp
                    break

                del tmp
            except Exception as exc:
                last_error = exc
                continue

        if np.isnan(pred_value) and last_error is not None:
            print(f"Warning: {sample_id} {canonical_clock} failed: {last_error}")

        results[canonical_clock] = pred_value

    del df, cpg_col, beta_vals, sample_matrix, adata
    gc.collect()
    return results


def infer_sample_id(bed_file: Path) -> str:
    name = bed_file.name
    if name.endswith(".dammet.bed"):
        return name[: -len(".dammet.bed")]
    if name.endswith(".bed"):
        return name[: -len(".bed")]
    return bed_file.stem


def default_illumina_path(bed_file: Path, sample_id: str) -> Path:
    return bed_file.with_name(f"{sample_id}_illumina.csv")


def default_predictions_path(bed_file: Path, sample_id: str) -> Path:
    return bed_file.with_name(f"{sample_id}_pyaging_predictions.csv")


def default_combined_predictions_path(first_bed_file: Path) -> Path:
    return first_bed_file.with_name("combined_pyaging_predictions.csv")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Convert one or more BED files to Illumina format and predict biological age with pyaging.",
    )
    parser.add_argument("bed_files", nargs="+", type=Path, help="One or more input BED files.")
    parser.add_argument(
        "--sample-id",
        help="Sample identifier to use for a single input BED file. Defaults to the BED filename.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=EPIC_MANIFEST,
        help=f"Path to the EPIC manifest CSV. Default: {EPIC_MANIFEST}",
    )
    parser.add_argument(
        "--illumina-out",
        type=Path,
        help="Output path for the generated Illumina CSV. Defaults to <sample_id>_illumina.csv next to the BED file.",
    )
    parser.add_argument(
        "--predictions-out",
        type=Path,
        help="Output path for the predictions CSV when one BED file is given. Defaults to <sample_id>_pyaging_predictions.csv next to the BED file.",
    )
    parser.add_argument(
        "--combined-out",
        type=Path,
        help="Output path for a combined predictions CSV when multiple BED files are given. Defaults to combined_pyaging_predictions.csv next to the first BED file.",
    )
    parser.add_argument(
        "--reuse-illumina",
        action="store_true",
        help="Reuse an existing Illumina CSV if present instead of rebuilding it.",
    )
    return parser


def process_bed_file(
    bed_file: Path,
    sample_id: str,
    manifest_path: Path,
    coord_to_cpg: Dict[str, str],
    reuse_illumina: bool,
    illumina_out: Path | None = None,
    predictions_out: Path | None = None,
) -> pd.DataFrame:
    resolved_bed_file = bed_file.resolve()
    if not resolved_bed_file.exists():
        raise FileNotFoundError(f"BED file not found: {resolved_bed_file}")

    resolved_illumina_out = (illumina_out or default_illumina_path(resolved_bed_file, sample_id)).resolve()
    resolved_predictions_out = (predictions_out or default_predictions_path(resolved_bed_file, sample_id)).resolve()

    resolved_illumina_out.parent.mkdir(parents=True, exist_ok=True)
    resolved_predictions_out.parent.mkdir(parents=True, exist_ok=True)

    if not reuse_illumina or not resolved_illumina_out.exists() or resolved_illumina_out.stat().st_size == 0:
        convert_bed_to_illumina_with_manifest(resolved_bed_file, resolved_illumina_out, sample_id, coord_to_cpg)

    predictions = predict_one_sample(resolved_illumina_out, sample_id)
    result_df = pd.DataFrame(
        [
            {
                "sample_id": sample_id,
                "bed_file": str(resolved_bed_file),
                "epic_manifest": str(manifest_path.resolve()),
                "illumina_file": str(resolved_illumina_out),
                "clock": clock,
                "predicted_age": predicted_age,
            }
            for clock, predicted_age in predictions.items()
        ]
    )
    result_df.to_csv(resolved_predictions_out, index=False)

    print(f"sample_id={sample_id}")
    print(f"bed_file={resolved_bed_file}")
    print(f"epic_manifest={manifest_path.resolve()}")
    print(f"illumina_file={resolved_illumina_out}")
    print(f"predictions_file={resolved_predictions_out}")
    print(result_df[["clock", "predicted_age"]].to_string(index=False))

    return result_df


def main() -> None:
    args = build_arg_parser().parse_args()
    bed_files = [bed_file.resolve() for bed_file in args.bed_files]

    if args.sample_id and len(bed_files) != 1:
        raise ValueError("--sample-id can only be used when exactly one BED file is provided.")

    if args.illumina_out and len(bed_files) != 1:
        raise ValueError("--illumina-out can only be used when exactly one BED file is provided.")

    if args.predictions_out and len(bed_files) != 1:
        raise ValueError("--predictions-out can only be used when exactly one BED file is provided.")

    ann = load_epic_manifest(args.manifest)
    coord_to_cpg = dict(zip(ann["coord_key"].values, ann["CpG"].values))

    all_results: List[pd.DataFrame] = []
    for index, bed_file in enumerate(bed_files):
        sample_id = args.sample_id if index == 0 and args.sample_id else infer_sample_id(bed_file)
        result_df = process_bed_file(
            bed_file=bed_file,
            sample_id=sample_id,
            manifest_path=args.manifest,
            coord_to_cpg=coord_to_cpg,
            reuse_illumina=args.reuse_illumina,
            illumina_out=args.illumina_out if len(bed_files) == 1 else None,
            predictions_out=args.predictions_out if len(bed_files) == 1 else None,
        )
        all_results.append(result_df)

    if len(all_results) > 1 or args.combined_out:
        combined_out = (args.combined_out or default_combined_predictions_path(bed_files[0])).resolve()
        combined_out.parent.mkdir(parents=True, exist_ok=True)
        combined_df = pd.concat(all_results, ignore_index=True)
        combined_df.to_csv(combined_out, index=False)
        print(f"combined_predictions_file={combined_out}")


if __name__ == "__main__":
    main()
