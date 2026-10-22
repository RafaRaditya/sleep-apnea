from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd
import numpy as np

from loaders import (
    DEFAULT_ROOT,
    _read_subject_details,
    info_to_row,
    inspect_record,
    list_records,
    load_record,
)


def build_inventory(root: Path) -> tuple[pd.DataFrame, dict]:
    subject_details = _read_subject_details(root)
    rows = []
    warnings = []

    for dataset in ["apnea-ecg", "slpdb", "ucddb"]:
        record_ids = list_records(dataset, root)
        print(f"[{dataset}] {len(record_ids)} records")
        for record_id in record_ids:
            info = inspect_record(dataset, record_id, root, subject_details)
            rows.append(info_to_row(info))

            if info.ecg_channel is None:
                warnings.append(f"{dataset}/{record_id}: no ECG channel")
            if info.resp_channel is None:
                warnings.append(f"{dataset}/{record_id}: no respiratory channel")
            if info.n_label_units == 0:
                warnings.append(f"{dataset}/{record_id}: no labels")
            if info.duration_s <= 0:
                warnings.append(f"{dataset}/{record_id}: non-positive duration")

    df = pd.DataFrame(rows)
    summary = summarize(df)
    summary["warnings"] = warnings
    return df, summary


def summarize(df: pd.DataFrame) -> dict:
    summary = {}
    for dataset, group in df.groupby("dataset", sort=False):
        n_events = pd.to_numeric(group["n_events"], errors="coerce")
        n_apnea_events = pd.to_numeric(group["n_apnea_events"], errors="coerce")
        summary[dataset] = {
            "n_records": int(len(group)),
            "n_subjects": int(group["subject_id"].nunique()),
            "total_duration_h": round(float(group["duration_s"].sum()) / 3600, 2),
            "n_label_units": int(group["n_label_units"].sum()),
            "n_apnea_units": int(group["n_apnea_units"].sum()),
            "apnea_prevalence": round(
                float(group["n_apnea_units"].sum()) / float(group["n_label_units"].sum()),
                4,
            ),
            "n_events": None if n_events.isna().all() else int(n_events.sum()),
            "n_apnea_events": None
            if n_apnea_events.isna().all()
            else int(n_apnea_events.sum()),
            "mean_ahi": None
            if group["ahi"].isna().all()
            else round(float(group["ahi"].mean()), 2),
            "mean_sleep_fraction": None
            if group["sleep_fraction"].isna().all()
            else round(float(group["sleep_fraction"].mean()), 4),
        }

    summary["overall"] = {
        "n_records": int(len(df)),
        "n_subjects": int(df["subject_id"].nunique()),
        "total_duration_h": round(float(df["duration_s"].sum()) / 3600, 2),
        "n_label_units": int(df["n_label_units"].sum()),
        "n_apnea_units": int(df["n_apnea_units"].sum()),
        "apnea_prevalence": round(
            float(df["n_apnea_units"].sum()) / float(df["n_label_units"].sum()), 4
        ),
    }
    return summary


def run_self_test(root: Path) -> int:
    failures = 0
    for dataset in ["apnea-ecg", "slpdb", "ucddb"]:
        record_id = list_records(dataset, root)[0]
        print(f"\n=== {dataset}/{record_id} ===")
        try:
            data = load_record(dataset, record_id, root)
        except Exception as exc:
            print(f"FAILED: {exc}")
            failures += 1
            continue

        print(f"subject={data.info.subject_id} duration_h={data.info.duration_s / 3600:.2f}")
        for name, signal in data.signals.items():
            fs = data.sampling_rates[name]
            print(
                f"  {name:4s} shape={signal.shape} fs={fs:.1f} "
                f"min={float(np.nanmin(signal)):.3f} max={float(np.nanmax(signal)):.3f} "
                f"nan={int(np.isnan(signal).sum())}"
            )
        print(
            f"  labels={len(data.labels)} positive={int(data.labels.sum())} "
            f"prevalence={data.labels.mean():.3f} unit={data.info.label_unit_s:.0f}s"
        )
        if data.sleep_stage is not None:
            print(f"  sleep_stage unique={sorted(set(map(str, data.sleep_stage)))[:8]}")
        if data.events is not None and not data.events.empty:
            print(data.events.head(3).to_string(index=False))
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Build dataset inventory for the sleep-apnea project.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    out_dir = args.out or root / "outputs"
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.self_test:
        failures = run_self_test(root)
        if failures:
            print(f"\n{failures} self-test(s) failed")
            return 1
        print("\nAll self-tests passed")
        return 0

    df, summary = build_inventory(root)

    csv_path = out_dir / "dataset_inventory.csv"
    json_path = out_dir / "dataset_summary.json"
    df.to_csv(csv_path, index=False, encoding="utf-8-sig")
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"\nWrote {csv_path}")
    print(f"Wrote {json_path}")
    print("\nSummary:")
    print(json.dumps({k: v for k, v in summary.items() if k != "warnings"}, indent=2))

    warnings = summary.get("warnings", [])
    if warnings:
        print(f"\nWarnings ({len(warnings)}):")
        for warning in warnings:
            print(f"  - {warning}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
