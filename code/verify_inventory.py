from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd
import pyedflib
import wfdb


APNEA_EVENTS = {"H", "HA", "OA", "X", "CA", "CAA"}
DATASETS = ["apnea-ecg", "slpdb", "ucddb"]
DATA_DIRS = {
    "apnea-ecg": "apnea-ecg/DESTINATION",
    "slpdb": "slpdb/DESTINATION",
    "ucddb": "ucddb/DESTINATION",
}


def parse_header(path: Path) -> dict:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    first = lines[0].split()
    fs = float(first[2].split("/")[0])
    n_samples = int(first[3])
    names = []
    for line in lines[1:]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split(maxsplit=8)
        if len(parts) >= 9:
            names.append(parts[8].strip())
    return {
        "n_signals": int(first[1]),
        "fs": fs,
        "n_samples": n_samples,
        "duration_s": n_samples / fs if fs else 0.0,
        "names": names,
    }


def clean_aux(aux: str) -> str:
    return aux.replace("\x00", "").strip().strip("'").strip()


def parse_ucddb_events(path: Path, start_dt: datetime, duration_s: float) -> list[tuple]:
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        tokens = line.split()
        if len(tokens) < 3 or not re.fullmatch(r"\d{2}:\d{2}:\d{2}", tokens[0]):
            continue
        event_type = tokens[1]
        j = 2
        if not tokens[j].isdigit():
            j += 1
        if j >= len(tokens) or not tokens[j].isdigit():
            continue
        event_duration = int(tokens[j])

        clock = datetime.strptime(tokens[0], "%H:%M:%S").time()
        event_dt = datetime.combine(start_dt.date(), clock)
        if event_dt < start_dt:
            event_dt += timedelta(days=1)
        onset = (event_dt - start_dt).total_seconds()
        if 0 <= onset <= duration_s + 60:
            rows.append(
                (
                    onset,
                    event_duration,
                    event_type,
                    event_type.upper().startswith(("HYP", "APNEA")),
                )
            )
    return rows


def mark_events(rows: list[tuple], n_epochs: int, epoch_s: float) -> np.ndarray:
    labels = np.zeros(n_epochs, dtype=np.int8)
    for onset, duration, _, is_apnea in rows:
        if not is_apnea:
            continue
        start = max(0, int(onset // epoch_s))
        end = min(n_epochs - 1, int((onset + duration - 1e-9) // epoch_s))
        labels[start : end + 1] = 1
    return labels


def checksum_map(data_dir: Path) -> dict[str, str]:
    path = data_dir / "SHA256SUMS.txt"
    if not path.exists():
        return {}
    result = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.split(maxsplit=1)
        if len(parts) == 2:
            result[parts[1].strip()] = parts[0]
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def verify_checksums(
    dataset: str,
    record_id: str,
    data_dir: Path,
    hashes: dict[str, str],
) -> list[str]:
    base = record_id[:-2] if dataset == "apnea-ecg" and record_id.endswith("er") else record_id
    if dataset == "apnea-ecg":
        names = [
            f"{record_id}.hea",
            f"{record_id}.apn",
            f"{base}.hea",
            f"{base}.dat",
            f"{base}r.hea",
            f"{base}r.dat",
        ]
    elif dataset == "slpdb":
        names = [f"{record_id}.hea", f"{record_id}.dat", f"{record_id}.st"]
    else:
        names = [
            f"{record_id}.rec",
            f"{record_id}_stage.txt",
            f"{record_id}_respevt.txt",
        ]

    problems = []
    for name in names:
        path = data_dir / name
        if not path.exists():
            problems.append(f"missing file {name}")
            continue
        expected = hashes.get(name)
        if expected is None:
            problems.append(f"no checksum entry for {name}")
            continue
        actual = sha256(path)
        if actual != expected:
            problems.append(f"checksum mismatch {name}")
    return problems


def compare(
    dataset: str,
    row: pd.Series,
    data_dir: Path,
) -> tuple[list[str], dict]:
    record_id = str(row["record_id"])
    problems: list[str] = []
    raw: dict = {}

    if dataset == "apnea-ecg":
        base = record_id[:-2] if record_id.endswith("er") else record_id
        header = parse_header(data_dir / f"{record_id}.hea")
        ann = wfdb.rdann(str(data_dir / record_id), "apn")
        n_label_units = len(ann.symbol)
        n_apnea_units = int(sum(1 for symbol in ann.symbol if symbol == "A"))

        raw.update(
            {
                "n_signals": header["n_signals"],
                "duration_s": header["duration_s"],
                "fs_ecg": header["fs"],
                "fs_resp": header["fs"],
                "n_label_units": n_label_units,
                "n_apnea_units": n_apnea_units,
                "n_events": None,
                "n_apnea_events": None,
            }
        )

    elif dataset == "slpdb":
        header = parse_header(data_dir / f"{record_id}.hea")
        ann = wfdb.rdann(str(data_dir / record_id), "st")
        n_label_units = len(ann.aux_note)
        n_apnea_units = 0
        n_events = 0
        n_apnea_events = 0
        for aux in ann.aux_note:
            tokens = clean_aux(aux).split()
            events = tokens[1:]
            n_events += len(events)
            has_apnea = False
            for token in events:
                if token in APNEA_EVENTS:
                    n_apnea_events += 1
                    has_apnea = True
            if has_apnea:
                n_apnea_units += 1

        raw.update(
            {
                "n_signals": header["n_signals"],
                "duration_s": header["duration_s"],
                "fs_ecg": header["fs"],
                "fs_resp": header["fs"],
                "n_label_units": n_label_units,
                "n_apnea_units": n_apnea_units,
                "n_events": n_events,
                "n_apnea_events": n_apnea_events,
            }
        )

    else:
        rec_path = data_dir / f"{record_id}.rec"
        reader = pyedflib.EdfReader(str(rec_path))
        try:
            signal_names = list(reader.getSignalLabels())
            fs_by_name = {
                name: float(reader.getSampleFrequency(i))
                for i, name in enumerate(signal_names)
            }
            duration_s = float(reader.getFileDuration())
            start_dt = reader.getStartdatetime()
        finally:
            reader.close()

        stage = [
            int(v)
            for v in (data_dir / f"{record_id}_stage.txt")
            .read_text(encoding="utf-8", errors="replace")
            .split()
        ]
        events = parse_ucddb_events(
            data_dir / f"{record_id}_respevt.txt", start_dt, duration_s
        )
        labels = mark_events(events, len(stage), 30.0)
        n_events = len(events)
        n_apnea_events = sum(1 for _, _, _, is_apnea in events if is_apnea)

        raw.update(
            {
                "n_signals": len(signal_names),
                "duration_s": duration_s,
                "fs_ecg": fs_by_name.get("ECG"),
                "fs_resp": fs_by_name.get("Flow"),
                "n_label_units": len(stage),
                "n_apnea_units": int(labels.sum()),
                "n_events": n_events,
                "n_apnea_events": n_apnea_events,
            }
        )

    exact_fields = ["n_signals", "n_label_units", "n_apnea_units"]
    for field in exact_fields:
        expected = row[field]
        if pd.isna(expected):
            continue
        if int(expected) != int(raw[field]):
            problems.append(f"{field}: csv={int(expected)} raw={int(raw[field])}")

    for field in ["n_events", "n_apnea_events"]:
        expected = row[field]
        if pd.isna(expected):
            continue
        if raw[field] is None or int(expected) != int(raw[field]):
            problems.append(f"{field}: csv={expected} raw={raw[field]}")

    for field, tol in [("duration_s", 0.5), ("fs_ecg", 0.01), ("fs_resp", 0.01)]:
        expected = row[field]
        if pd.isna(expected):
            continue
        actual = raw[field]
        if actual is None or abs(float(expected) - float(actual)) > tol:
            problems.append(f"{field}: csv={expected} raw={actual}")

    expected_prev = row["apnea_prevalence"]
    if not pd.isna(expected_prev):
        actual_prev = raw["n_apnea_units"] / raw["n_label_units"] if raw["n_label_units"] else 0.0
        if abs(float(expected_prev) - actual_prev) > 0.001:
            problems.append(
                f"apnea_prevalence: csv={expected_prev} raw={actual_prev:.4f}"
            )

    if not row["ecg_channel"] or pd.isna(row["ecg_channel"]):
        problems.append("ecg_channel missing")
    if not row["resp_channel"] or pd.isna(row["resp_channel"]):
        problems.append("resp_channel missing")

    return problems, raw


def verify_json(csv: pd.DataFrame, summary: dict) -> list[str]:
    problems = []
    for dataset in DATASETS:
        group = csv[csv["dataset"] == dataset]
        expected = summary.get(dataset, {})
        checks = {
            "n_records": len(group),
            "n_subjects": group["subject_id"].nunique(),
            "n_label_units": int(group["n_label_units"].sum()),
            "n_apnea_units": int(group["n_apnea_units"].sum()),
        }
        for key, value in checks.items():
            if int(expected.get(key, -1)) != int(value):
                problems.append(f"json {dataset}.{key}: json={expected.get(key)} csv={value}")

        expected_prev = float(group["n_apnea_units"].sum()) / float(
            group["n_label_units"].sum()
        )
        if abs(float(expected.get("apnea_prevalence", -1)) - expected_prev) > 0.001:
            problems.append(
                f"json {dataset}.apnea_prevalence: json={expected.get('apnea_prevalence')} csv={expected_prev:.4f}"
            )

    overall = summary.get("overall", {})
    overall_checks = {
        "n_records": len(csv),
        "n_subjects": csv["subject_id"].nunique(),
        "n_label_units": int(csv["n_label_units"].sum()),
        "n_apnea_units": int(csv["n_apnea_units"].sum()),
    }
    for key, value in overall_checks.items():
        if int(overall.get(key, -1)) != int(value):
            problems.append(f"json overall.{key}: json={overall.get(key)} csv={value}")
    overall_prev = float(csv["n_apnea_units"].sum()) / float(csv["n_label_units"].sum())
    if abs(float(overall.get("apnea_prevalence", -1)) - overall_prev) > 0.001:
        problems.append(
            f"json overall.apnea_prevalence: json={overall.get('apnea_prevalence')} csv={overall_prev:.4f}"
        )
    return problems


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Independently verify dataset_inventory.csv and dataset_summary.json."
    )
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--checksums",
        choices=["sample", "all", "none"],
        default="sample",
        help="Verify raw-file SHA256 against official SHA256SUMS.txt.",
    )
    args = parser.parse_args()

    root = args.root.resolve()
    out_dir = root / "outputs"
    csv = pd.read_csv(out_dir / "dataset_inventory.csv")
    summary = json.loads((out_dir / "dataset_summary.json").read_text(encoding="utf-8"))

    results = []
    all_problems = []

    print("Independent recomputation from raw files...")
    for _, row in csv.iterrows():
        dataset = str(row["dataset"])
        record_id = str(row["record_id"])
        data_dir = root / DATA_DIRS[dataset]
        problems, raw = compare(dataset, row, data_dir)
        results.append(
            {
                "dataset": dataset,
                "record_id": record_id,
                "status": "PASS" if not problems else "FAIL",
                "mismatches": "; ".join(problems),
                "raw_n_label_units": raw.get("n_label_units"),
                "raw_n_apnea_units": raw.get("n_apnea_units"),
            }
        )
        if problems:
            all_problems.append(f"{dataset}/{record_id}: " + "; ".join(problems))

    print("Checking JSON aggregates against CSV...")
    json_problems = verify_json(csv, summary)
    all_problems.extend(json_problems)

    checksum_problems = []
    if args.checksums != "none":
        if args.checksums == "sample":
            sample_ids = csv.groupby("dataset", sort=False).head(1)
            check_rows = sample_ids
        else:
            check_rows = csv

        print(f"Verifying raw-file checksums ({args.checksums})...")
        for _, row in check_rows.iterrows():
            dataset = str(row["dataset"])
            record_id = str(row["record_id"])
            data_dir = root / DATA_DIRS[dataset]
            hashes = checksum_map(data_dir)
            problems = verify_checksums(dataset, record_id, data_dir, hashes)
            checksum_problems.extend(
                [f"{dataset}/{record_id}: {problem}" for problem in problems]
            )

    result_df = pd.DataFrame(results)
    result_df.to_csv(out_dir / "verification_results.csv", index=False, encoding="utf-8-sig")

    n_pass = int((result_df["status"] == "PASS").sum())
    n_fail = int((result_df["status"] == "FAIL").sum())

    report = []
    report.append("# Verification Report — Inventory vs Raw Files\n")
    report.append("## Ringkasan\n")
    report.append(f"- Record diverifikasi ulang: {len(result_df)}")
    report.append(f"- PASS: {n_pass}")
    report.append(f"- FAIL: {n_fail}")
    report.append(f"- Masalah agregat JSON: {len(json_problems)}")
    report.append(f"- Masalah checksum ({args.checksums}): {len(checksum_problems)}")
    report.append("")
    if all_problems or checksum_problems:
        report.append("## Mismatch\n")
        for problem in all_problems + checksum_problems:
            report.append(f"- {problem}")
    else:
        report.append("## Hasil\n")
        report.append("Semua nilai CSV/JSON cocok dengan hasil hitung ulang dari file mentah, dan checksum file mentah cocok dengan `SHA256SUMS.txt` resmi.")
    (out_dir / "verification_report.md").write_text("\n".join(report), encoding="utf-8")

    print(f"\nRecords: {len(result_df)} | PASS: {n_pass} | FAIL: {n_fail}")
    print(f"JSON aggregate problems: {len(json_problems)}")
    print(f"Checksum problems ({args.checksums}): {len(checksum_problems)}")
    print("Wrote outputs/verification_results.csv and outputs/verification_report.md")

    return 0 if not (all_problems or checksum_problems) else 1


if __name__ == "__main__":
    raise SystemExit(main())
