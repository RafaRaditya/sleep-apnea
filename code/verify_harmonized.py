from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


DATASETS = ["apnea-ecg", "slpdb", "ucddb"]


def check_record(row: pd.Series, harm_dir: Path) -> list[str]:
    problems = []
    stem = f"{row['dataset']}__{row['record_id']}"
    path = harm_dir / "signals" / f"{stem}.npz"
    if not path.exists():
        return [f"missing {path.name}"]

    data = np.load(path, allow_pickle=False)
    required = [
        "ecg",
        "resp",
        "fs",
        "labels_30s",
        "label_mask_30s",
        "label_times_30s",
        "sleep_stage_30s",
        "native_labels",
        "native_label_times_s",
        "native_label_unit_s",
        "meta_json",
    ]
    for key in required:
        if key not in data.files:
            problems.append(f"missing key {key}")
    if problems:
        return problems

    ecg = data["ecg"]
    resp = data["resp"]
    fs = float(data["fs"])
    labels = data["labels_30s"]
    mask = data["label_mask_30s"]
    sleep_stage = data["sleep_stage_30s"]
    n_samples = int(row["n_samples"])
    n_epochs = int(row["n_epochs_30s"])

    if ecg.ndim != 1 or resp.ndim != 1:
        problems.append("ecg/resp not 1-D")
    if len(ecg) != n_samples or len(resp) != n_samples:
        problems.append(
            f"length: ecg={len(ecg)} resp={len(resp)} manifest={n_samples}"
        )
    if not np.isfinite(ecg).all():
        problems.append("ecg contains NaN/inf")
    if not np.isfinite(resp).all():
        problems.append("resp contains NaN/inf")
    if abs(fs - float(row["target_fs"])) > 1e-3:
        problems.append(f"fs: npz={fs} manifest={row['target_fs']}")

    lengths = {
        "labels_30s": len(labels),
        "label_mask_30s": len(mask),
        "label_times_30s": len(data["label_times_30s"]),
        "sleep_stage_30s": len(sleep_stage),
    }
    for name, length in lengths.items():
        if length != n_epochs:
            problems.append(f"{name} length={length} manifest={n_epochs}")

    expected_epochs = n_samples // int(round(30 * fs))
    if n_epochs != expected_epochs:
        problems.append(f"n_epochs={n_epochs} expected={expected_epochs}")
    if int(mask.sum()) != int(row["n_labeled_30s"]):
        problems.append(
            f"label_mask sum={int(mask.sum())} manifest={int(row['n_labeled_30s'])}"
        )
    if int((labels * mask).sum()) != int(row["n_apnea_30s"]):
        problems.append(
            f"positive labels={int((labels * mask).sum())} manifest={int(row['n_apnea_30s'])}"
        )
    if not set(np.unique(labels)).issubset({0, 1}):
        problems.append(f"labels values={np.unique(labels).tolist()}")
    if not set(np.unique(mask)).issubset({0, 1}):
        problems.append(f"mask values={np.unique(mask).tolist()}")
    if not set(np.unique(sleep_stage)).issubset(set(range(-1, 7))):
        problems.append(f"sleep_stage values={np.unique(sleep_stage).tolist()}")

    duration_h = n_samples / fs / 3600
    if abs(duration_h - float(row["duration_h"])) > 0.01:
        problems.append(
            f"duration_h={duration_h:.4f} manifest={row['duration_h']}"
        )

    native_labels = data["native_labels"]
    if len(native_labels) != int(row["n_native_label_units"]):
        problems.append(
            f"native_labels length={len(native_labels)} manifest={int(row['n_native_label_units'])}"
        )
    if abs(float(data["native_label_unit_s"]) - float(row["native_label_unit_s"])) > 1e-3:
        problems.append("native_label_unit_s mismatch")

    if row["dataset"] == "ucddb":
        n_events = len(data["event_onset_s"])
        if len(data["event_duration_s"]) != n_events:
            problems.append("event_duration_s length mismatch")
        if len(data["event_apnea"]) != n_events:
            problems.append("event_apnea length mismatch")
        if len(data["event_type"]) != n_events:
            problems.append("event_type length mismatch")
    elif row["dataset"] == "slpdb":
        if len(data["slpdb_event_tokens"]) != int(row["n_native_label_units"]):
            problems.append("slpdb_event_tokens length mismatch")
        if len(data["slpdb_event_epoch"]) != int(row["n_native_label_units"]):
            problems.append("slpdb_event_epoch length mismatch")

    try:
        meta = json.loads(str(data["meta_json"].item()))
        if int(meta["harmonized"]["n_samples"]) != n_samples:
            problems.append("meta n_samples mismatch")
        if meta["source"]["resp_channel"] != row["resp_channel"]:
            problems.append("meta resp_channel mismatch")
    except Exception as exc:
        problems.append(f"meta_json parse failed: {exc}")

    return problems


def verify_summary(manifest: pd.DataFrame, summary: dict) -> list[str]:
    problems = []
    for dataset in DATASETS:
        group = manifest[manifest["dataset"] == dataset]
        if group.empty:
            continue
        expected = summary.get(dataset, {})
        n_labeled = int(group["n_labeled_30s"].sum())
        n_apnea = int(group["n_apnea_30s"].sum())
        checks = {
            "n_records": len(group),
            "n_subjects": group["subject_id"].nunique(),
            "n_epochs_30s": int(group["n_epochs_30s"].sum()),
            "n_labeled_30s": n_labeled,
            "n_apnea_30s": n_apnea,
        }
        for key, value in checks.items():
            if int(expected.get(key, -1)) != int(value):
                problems.append(f"summary {dataset}.{key}: {expected.get(key)} != {value}")
        prevalence = n_apnea / n_labeled if n_labeled else 0.0
        if abs(float(expected.get("apnea_prevalence_30s", -1)) - prevalence) > 0.001:
            problems.append(
                f"summary {dataset}.apnea_prevalence_30s: "
                f"{expected.get('apnea_prevalence_30s')} != {prevalence:.4f}"
            )
        total_h = float(group["duration_h"].sum())
        if abs(float(expected.get("total_duration_h", -1)) - total_h) > 0.02:
            problems.append(
                f"summary {dataset}.total_duration_h: "
                f"{expected.get('total_duration_h')} != {total_h:.2f}"
            )

    overall = summary.get("overall", {})
    n_labeled = int(manifest["n_labeled_30s"].sum())
    n_apnea = int(manifest["n_apnea_30s"].sum())
    overall_checks = {
        "n_records": len(manifest),
        "n_subjects": manifest["subject_id"].nunique(),
        "n_epochs_30s": int(manifest["n_epochs_30s"].sum()),
        "n_labeled_30s": n_labeled,
        "n_apnea_30s": n_apnea,
    }
    for key, value in overall_checks.items():
        if int(overall.get(key, -1)) != int(value):
            problems.append(f"summary overall.{key}: {overall.get(key)} != {value}")
    prevalence = n_apnea / n_labeled if n_labeled else 0.0
    if abs(float(overall.get("apnea_prevalence_30s", -1)) - prevalence) > 0.001:
        problems.append(
            f"summary overall.apnea_prevalence_30s: "
            f"{overall.get('apnea_prevalence_30s')} != {prevalence:.4f}"
        )
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify harmonized .npz files and manifest.")
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--harmonized", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    harm_dir = args.harmonized or root / "outputs" / "harmonized"
    manifest = pd.read_csv(harm_dir / "manifest.csv")
    summary = json.loads((harm_dir / "summary.json").read_text(encoding="utf-8"))

    results = []
    problems_all = []
    for _, row in manifest.iterrows():
        problems = check_record(row, harm_dir)
        results.append(
            {
                "dataset": row["dataset"],
                "record_id": row["record_id"],
                "status": "PASS" if not problems else "FAIL",
                "mismatches": "; ".join(problems),
            }
        )
        if problems:
            problems_all.append(
                f"{row['dataset']}/{row['record_id']}: " + "; ".join(problems)
            )

    summary_problems = verify_summary(manifest, summary)
    problems_all.extend(summary_problems)

    result_df = pd.DataFrame(results)
    result_df.to_csv(
        harm_dir / "verification_results.csv", index=False, encoding="utf-8-sig"
    )

    n_pass = int((result_df["status"] == "PASS").sum())
    n_fail = int((result_df["status"] == "FAIL").sum())

    report = [
        "# Verification Report — Harmonized Dataset\n",
        "## Ringkasan\n",
        f"- Record diverifikasi: {len(result_df)}",
        f"- PASS: {n_pass}",
        f"- FAIL: {n_fail}",
        f"- Masalah agregat summary: {len(summary_problems)}",
        "",
    ]
    if problems_all:
        report.append("## Mismatch\n")
        for problem in problems_all:
            report.append(f"- {problem}")
    else:
        report.append("## Hasil\n")
        report.append(
            "Semua `.npz` konsisten dengan `manifest.csv`, dan `summary.json` konsisten dengan manifest."
        )
    (harm_dir / "verification_report.md").write_text(
        "\n".join(report), encoding="utf-8"
    )

    print(f"Records: {len(result_df)} | PASS: {n_pass} | FAIL: {n_fail}")
    print(f"Summary problems: {len(summary_problems)}")
    print("Wrote outputs/harmonized/verification_report.md")
    return 0 if not problems_all else 1


if __name__ == "__main__":
    raise SystemExit(main())
