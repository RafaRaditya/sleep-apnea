from __future__ import annotations

import argparse
import json
import math
from fractions import Fraction
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.signal import resample_poly

from loaders import DEFAULT_ROOT, list_records, load_record


DATASETS = ["apnea-ecg", "slpdb", "ucddb"]
EPOCH_S = 30.0

RESP_PRIORITY = [
    "Resp N",
    "Resp (nasal)",
    "Flow",
    "Resp (sum)",
    "Resp C",
    "Resp (chest)",
    "Resp A",
    "Resp (abdominal)",
    "Resp (abdomen)",
]

SLPDB_STAGE_MAP = {
    "W": 0,
    "1": 1,
    "2": 2,
    "3": 3,
    "4": 4,
    "R": 5,
    "MT": 6,
    "M": 6,
}
UCDDB_STAGE_MAP = {
    0: 0,  # Wake
    1: 5,  # REM
    2: 1,  # Stage 1
    3: 2,  # Stage 2
    4: 3,  # Stage 3
    5: 4,  # Stage 4
    6: 6,  # Artifact
    7: -1,  # Indeterminate
}


def max_true_run(mask: np.ndarray) -> int:
    if not mask.any():
        return 0
    padded = np.concatenate(([False], mask, [False]))
    diff = np.diff(padded.astype(np.int8))
    starts = np.where(diff == 1)[0]
    ends = np.where(diff == -1)[0]
    return int((ends - starts).max()) if len(starts) else 0


def interpolate_nans(x: np.ndarray, fs: float) -> tuple[np.ndarray, float, float]:
    x = np.asarray(x, dtype=np.float64)
    mask = ~np.isfinite(x)
    if not mask.any():
        return x.astype(np.float32), 0.0, 0.0

    nan_fraction = float(mask.mean())
    max_gap_s = max_true_run(mask) / fs if fs else 0.0
    valid = ~mask
    if not valid.any():
        return np.zeros_like(x, dtype=np.float32), 1.0, len(x) / fs if fs else float("inf")

    idx = np.arange(len(x))
    filled = np.interp(idx, idx[valid], x[valid])
    return filled.astype(np.float32), nan_fraction, max_gap_s


def resample_to(x: np.ndarray, fs_in: float, fs_out: float) -> np.ndarray:
    if abs(fs_in - fs_out) < 1e-9:
        return np.asarray(x, dtype=np.float32)
    ratio = Fraction(fs_out / fs_in).limit_denominator(1000)
    return resample_poly(np.asarray(x, dtype=np.float64), ratio.numerator, ratio.denominator).astype(
        np.float32
    )


def select_respiratory_channel(data, info) -> tuple[np.ndarray, str]:
    resp = np.asarray(data.signals["resp"])
    if resp.ndim == 1:
        return resp, info.resp_channel or "resp"

    names = (info.resp_channel or "").split(";")
    for target in RESP_PRIORITY:
        if target in names:
            return resp[:, names.index(target)], target
    return resp[:, 0], names[0] if names else "resp"


def map_labels_to_epochs(
    label_times_s: np.ndarray,
    labels: np.ndarray,
    unit_s: float,
    n_epochs: int,
) -> tuple[np.ndarray, np.ndarray]:
    labels_30 = np.zeros(n_epochs, dtype=np.int8)
    mask_30 = np.zeros(n_epochs, dtype=np.int8)
    for time_s, label in zip(label_times_s, labels):
        start = int(round(float(time_s) / EPOCH_S))
        end = int(round((float(time_s) + float(unit_s)) / EPOCH_S)) - 1
        if end < start:
            end = start
        start = max(0, start)
        end = min(n_epochs - 1, end)
        if start > end:
            continue
        labels_30[start : end + 1] = int(label)
        mask_30[start : end + 1] = 1
    return labels_30, mask_30


def harmonize_sleep_stage(stage, dataset: str, n_epochs: int) -> np.ndarray:
    out = np.full(n_epochs, -1, dtype=np.int8)
    if stage is None:
        return out

    m = min(len(stage), n_epochs)
    if dataset == "slpdb":
        for i in range(m):
            out[i] = SLPDB_STAGE_MAP.get(str(stage[i]).strip(), -1)
    elif dataset == "ucddb":
        for i in range(m):
            out[i] = UCDDB_STAGE_MAP.get(int(stage[i]), -1)
    return out


def harmonize_record(
    dataset: str,
    record_id: str,
    root: Path,
    target_fs: float,
) -> tuple[dict, dict, dict]:
    data = load_record(dataset, record_id, root)
    info = data.info

    ecg_raw = np.asarray(data.signals["ecg"], dtype=np.float32)
    resp_raw, resp_name = select_respiratory_channel(data, info)

    fs_ecg = float(data.sampling_rates["ecg"])
    fs_resp = float(data.sampling_rates["resp"])

    ecg_filled, nan_ecg, max_gap_ecg = interpolate_nans(ecg_raw, fs_ecg)
    resp_filled, nan_resp, max_gap_resp = interpolate_nans(resp_raw, fs_resp)

    ecg = resample_to(ecg_filled, fs_ecg, target_fs)
    resp = resample_to(resp_filled, fs_resp, target_fs)

    n = min(len(ecg), len(resp))
    ecg = ecg[:n]
    resp = resp[:n]

    n_epochs = int(math.floor(n / (EPOCH_S * target_fs)))
    if n_epochs <= 0:
        raise ValueError(f"{dataset}/{record_id}: no full 30-s epoch after resampling")

    native_labels = np.asarray(data.labels, dtype=np.int8)
    native_times = np.asarray(data.label_times_s, dtype=np.float64)
    labels_30, mask_30 = map_labels_to_epochs(
        native_times, native_labels, float(info.label_unit_s), n_epochs
    )
    sleep_stage = harmonize_sleep_stage(data.sleep_stage, dataset, n_epochs)

    event_onset = np.array([], dtype=np.float32)
    event_duration = np.array([], dtype=np.float32)
    event_apnea = np.array([], dtype=np.int8)
    event_type = np.array([], dtype="U16")
    slpdb_event_epoch = np.array([], dtype=np.int16)
    slpdb_event_tokens = np.array([], dtype="U32")

    if dataset == "ucddb" and data.events is not None and not data.events.empty:
        event_onset = data.events["onset_s"].to_numpy(dtype=np.float32)
        event_duration = data.events["duration_s"].to_numpy(dtype=np.float32)
        event_apnea = data.events["apnea"].to_numpy(dtype=np.int8)
        event_type = np.asarray(
            data.events["type"].astype(str).tolist(), dtype="U16"
        )
    elif dataset == "slpdb" and data.events is not None and not data.events.empty:
        slpdb_event_epoch = data.events["epoch"].to_numpy(dtype=np.int16)
        slpdb_event_tokens = np.asarray(
            data.events["events"].astype(str).tolist(), dtype="U32"
        )

    n_labeled = int(mask_30.sum())
    n_apnea = int((labels_30 * mask_30).sum())
    prevalence = n_apnea / n_labeled if n_labeled else 0.0
    duration_h = n / target_fs / 3600
    sleep_known = float((sleep_stage != -1).mean()) if n_epochs else 0.0
    if data.sleep_stage is None:
        n_sleep = None
        sleep_fraction = None
    else:
        is_sleep = np.isin(sleep_stage, [1, 2, 3, 4, 5])
        n_sleep = int(is_sleep.sum())
        sleep_fraction = float(is_sleep.mean()) if n_epochs else None

    meta = {
        "dataset": dataset,
        "record_id": record_id,
        "subject_id": info.subject_id,
        "source": {
            "label_source": info.label_source,
            "label_unit_s": info.label_unit_s,
            "signal_names": info.signal_names,
            "resp_channel": resp_name,
            "resp_channels_available": info.resp_channel,
            "native_fs_ecg": fs_ecg,
            "native_fs_resp": fs_resp,
        },
        "harmonized": {
            "target_fs": target_fs,
            "n_samples": int(n),
            "duration_h": round(duration_h, 4),
            "n_epochs_30s": int(n_epochs),
            "n_labeled_30s": n_labeled,
            "n_apnea_30s": n_apnea,
            "apnea_prevalence_30s": round(prevalence, 6),
            "nan_fraction_ecg": nan_ecg,
            "nan_fraction_resp": nan_resp,
            "max_nan_gap_ecg_s": max_gap_ecg,
            "max_nan_gap_resp_s": max_gap_resp,
            "sleep_known_fraction": round(sleep_known, 6),
            "n_sleep_30s": n_sleep,
            "sleep_fraction_30s": None
            if sleep_fraction is None
            else round(sleep_fraction, 6),
        },
        "notes": info.notes,
    }

    payload = {
        "ecg": ecg,
        "resp": resp,
        "fs": np.float32(target_fs),
        "labels_30s": labels_30,
        "label_mask_30s": mask_30,
        "label_times_30s": (np.arange(n_epochs) * EPOCH_S).astype(np.float32),
        "sleep_stage_30s": sleep_stage,
        "native_labels": native_labels,
        "native_label_times_s": native_times.astype(np.float32),
        "native_label_unit_s": np.float32(info.label_unit_s),
        "event_onset_s": event_onset,
        "event_duration_s": event_duration,
        "event_apnea": event_apnea,
        "event_type": event_type,
        "slpdb_event_epoch": slpdb_event_epoch,
        "slpdb_event_tokens": slpdb_event_tokens,
    }

    manifest_row = {
        "dataset": dataset,
        "record_id": record_id,
        "subject_id": info.subject_id,
        "native_fs_ecg": fs_ecg,
        "native_fs_resp": fs_resp,
        "target_fs": target_fs,
        "resp_channel": resp_name,
        "duration_h": round(duration_h, 3),
        "n_samples": int(n),
        "n_epochs_30s": int(n_epochs),
        "n_labeled_30s": n_labeled,
        "n_apnea_30s": n_apnea,
        "apnea_prevalence_30s": round(prevalence, 6),
        "n_native_label_units": int(len(native_labels)),
        "native_label_unit_s": float(info.label_unit_s),
        "nan_fraction_ecg": round(nan_ecg, 6),
        "nan_fraction_resp": round(nan_resp, 6),
        "max_nan_gap_resp_s": round(max_gap_resp, 2),
        "sleep_known_fraction": round(sleep_known, 6),
        "n_sleep_30s": n_sleep,
        "sleep_fraction_30s": None
        if sleep_fraction is None
        else round(sleep_fraction, 6),
        "notes": info.notes,
    }

    return meta, manifest_row, payload


def summarize(manifest: pd.DataFrame) -> dict:
    summary = {}
    for dataset, group in manifest.groupby("dataset", sort=False):
        n_labeled = int(group["n_labeled_30s"].sum())
        n_apnea = int(group["n_apnea_30s"].sum())
        summary[dataset] = {
            "n_records": int(len(group)),
            "n_subjects": int(group["subject_id"].nunique()),
            "total_duration_h": round(float(group["duration_h"].sum()), 2),
            "n_epochs_30s": int(group["n_epochs_30s"].sum()),
            "n_labeled_30s": n_labeled,
            "n_apnea_30s": n_apnea,
            "apnea_prevalence_30s": round(n_apnea / n_labeled, 4) if n_labeled else 0.0,
            "mean_nan_fraction_resp": round(float(group["nan_fraction_resp"].mean()), 6),
            "mean_sleep_known_fraction": round(
                float(group["sleep_known_fraction"].mean()), 4
            ),
            "mean_sleep_fraction_30s": None
            if group["sleep_fraction_30s"].isna().all()
            else round(float(group["sleep_fraction_30s"].mean()), 4),
        }

    n_labeled = int(manifest["n_labeled_30s"].sum())
    n_apnea = int(manifest["n_apnea_30s"].sum())
    summary["overall"] = {
        "n_records": int(len(manifest)),
        "n_subjects": int(manifest["subject_id"].nunique()),
        "total_duration_h": round(float(manifest["duration_h"].sum()), 2),
        "n_epochs_30s": int(manifest["n_epochs_30s"].sum()),
        "n_labeled_30s": n_labeled,
        "n_apnea_30s": n_apnea,
        "apnea_prevalence_30s": round(n_apnea / n_labeled, 4) if n_labeled else 0.0,
    }
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description="Harmonize the three sleep-apnea datasets.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--out", type=Path, default=None)
    parser.add_argument("--fs", type=float, default=100.0)
    parser.add_argument(
        "--datasets",
        type=str,
        default=",".join(DATASETS),
        help="Comma-separated subset, e.g. apnea-ecg,ucddb",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Max records per dataset (for a quick test).",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    root = args.root.resolve()
    out_dir = args.out or root / "outputs" / "harmonized"
    signal_dir = out_dir / "signals"
    signal_dir.mkdir(parents=True, exist_ok=True)

    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]
    manifest_rows = []
    failures = []

    total = 0
    for dataset in datasets:
        record_ids = list_records(dataset, root)
        if args.limit:
            record_ids = record_ids[: args.limit]
        total += len(record_ids)

    index = 0
    for dataset in datasets:
        record_ids = list_records(dataset, root)
        if args.limit:
            record_ids = record_ids[: args.limit]

        for record_id in record_ids:
            index += 1
            stem = f"{dataset}__{record_id}"
            npz_path = signal_dir / f"{stem}.npz"
            meta_path = signal_dir / f"{stem}.json"
            prefix = f"[{index}/{total}] {dataset}/{record_id}"

            if npz_path.exists() and meta_path.exists() and not args.overwrite:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                manifest_rows.append(meta["manifest"])
                print(f"{prefix}: skip (exists)")
                continue

            try:
                meta, manifest_row, payload = harmonize_record(
                    dataset, record_id, root, args.fs
                )
                payload["meta_json"] = np.array(json.dumps(meta, indent=2))
                np.savez_compressed(npz_path, **payload)
                meta["manifest"] = manifest_row
                meta_path.write_text(json.dumps(meta, indent=2), encoding="utf-8")
                manifest_rows.append(manifest_row)
                print(
                    f"{prefix}: {manifest_row['duration_h']:.2f} h | "
                    f"{manifest_row['n_epochs_30s']} epochs | "
                    f"apnea {manifest_row['apnea_prevalence_30s']:.3f}"
                )
            except Exception as exc:
                failures.append(f"{dataset}/{record_id}: {exc}")
                print(f"{prefix}: FAILED {exc}")

    manifest = pd.DataFrame(manifest_rows)
    manifest_path = out_dir / "manifest.csv"
    manifest.to_csv(manifest_path, index=False, encoding="utf-8-sig")

    if manifest.empty:
        summary = {
            "overall": {
                "n_records": 0,
                "n_subjects": 0,
                "total_duration_h": 0.0,
                "n_epochs_30s": 0,
                "n_labeled_30s": 0,
                "n_apnea_30s": 0,
                "apnea_prevalence_30s": 0.0,
            }
        }
    else:
        summary = summarize(manifest)
    summary["failures"] = failures
    (out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    print(f"\nWrote {manifest_path}")
    print(f"Wrote {out_dir / 'summary.json'}")
    print(json.dumps({k: v for k, v in summary.items() if k != "failures"}, indent=2))

    if failures:
        print(f"\nFailures ({len(failures)}):")
        for failure in failures:
            print(f"  - {failure}")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
