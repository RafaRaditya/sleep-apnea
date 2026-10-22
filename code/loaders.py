from __future__ import annotations

import json
import os
import re
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

try:
    import wfdb
except ImportError as exc:
    raise ImportError("Install wfdb: pip install wfdb") from exc

try:
    import pyedflib
except ImportError as exc:
    raise ImportError("Install pyedflib: pip install pyedflib") from exc


DEFAULT_ROOT = Path(
    os.environ.get("RESEARCH_ROOT", Path(__file__).resolve().parents[1])
)

DATA_DIRS = {
    "apnea-ecg": "apnea-ecg/DESTINATION",
    "slpdb": "slpdb/DESTINATION",
    "ucddb": "ucddb/DESTINATION",
}

APNEA_EVENTS = {"H", "HA", "OA", "X", "CA", "CAA"}
SLEEP_STAGES = {"1", "2", "3", "4", "R"}
SLPDB_RESP_PREFERENCE = [
    "Resp (sum)",
    "Resp (nasal)",
    "Resp (chest)",
    "Resp (abdominal)",
    "Resp (abdomen)",
]
UCDDB_RESP_PREFERENCE = ["Flow", "Sum", "ribcage", "abdo"]


@dataclass
class RecordInfo:
    dataset: str
    record_id: str
    subject_id: str
    n_signals: int
    signal_names: list[str]
    ecg_channel: str | None
    resp_channel: str | None
    fs_ecg: float | None
    fs_resp: float | None
    duration_s: float
    label_source: str
    label_unit_s: float
    n_label_units: int
    n_apnea_units: int
    n_events: int | None
    n_apnea_events: int | None
    event_types: dict[str, int] = field(default_factory=dict)
    ahi: float | None = None
    sleep_efficiency: float | None = None
    sleep_fraction: float | None = None
    notes: str = ""

    @property
    def apnea_prevalence(self) -> float:
        return self.n_apnea_units / self.n_label_units if self.n_label_units else 0.0


@dataclass
class RecordData:
    info: RecordInfo
    signals: dict[str, np.ndarray]
    sampling_rates: dict[str, float]
    label_times_s: np.ndarray
    labels: np.ndarray
    sleep_stage: np.ndarray | None = None
    events: pd.DataFrame | None = None


def _data_dir(dataset: str, root: Path | None = None) -> Path:
    root = Path(root) if root else DEFAULT_ROOT
    return root / DATA_DIRS[dataset]


def parse_wfdb_header(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    if not lines:
        raise ValueError(f"Empty header: {path}")

    first = lines[0].split()
    if len(first) < 4:
        raise ValueError(f"Malformed header: {path}")

    record_name = first[0]
    n_signals = int(first[1])
    fs = float(first[2].split("/")[0])
    n_samples = int(first[3])

    signals = []
    for line in lines[1:]:
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        parts = line.split(maxsplit=8)
        if len(parts) < 9:
            continue
        signals.append(
            {
                "file": parts[0],
                "format": parts[1],
                "description": parts[8].strip(),
            }
        )

    return {
        "record_name": record_name,
        "n_signals": n_signals,
        "fs": fs,
        "n_samples": n_samples,
        "duration_s": n_samples / fs if fs else 0.0,
        "signals": signals,
        "signal_names": [s["description"] for s in signals],
    }


def _clean_aux(aux: str) -> str:
    return aux.replace("\x00", "").strip().strip("'").strip()


def _check_spacing(samples: np.ndarray, fs: float, expected_s: float, name: str) -> list[str]:
    warnings = []
    if len(samples) > 1:
        spacing = np.diff(samples) / fs
        median = float(np.median(spacing))
        if abs(median - expected_s) > 0.5:
            warnings.append(
                f"{name}: label spacing {median:.2f}s != expected {expected_s:.0f}s"
            )
    return warnings


def _read_ucddb_stage(path: Path) -> np.ndarray:
    values = [int(v) for v in path.read_text(encoding="utf-8", errors="replace").split()]
    return np.asarray(values, dtype=np.int16)


def parse_ucddb_events(
    path: str | Path,
    start_datetime: datetime,
    duration_s: float,
) -> pd.DataFrame:
    path = Path(path)
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        tokens = line.split()
        if len(tokens) < 3:
            continue
        if not re.fullmatch(r"\d{2}:\d{2}:\d{2}", tokens[0]):
            continue

        event_type = tokens[1]
        duration_idx = 2
        if not tokens[duration_idx].isdigit():
            duration_idx += 1
        if duration_idx >= len(tokens) or not tokens[duration_idx].isdigit():
            continue
        event_duration = int(tokens[duration_idx])

        clock = datetime.strptime(tokens[0], "%H:%M:%S").time()
        event_dt = datetime.combine(start_datetime.date(), clock)
        if event_dt < start_datetime:
            event_dt += timedelta(days=1)
        onset_s = (event_dt - start_datetime).total_seconds()

        if onset_s < 0 or onset_s > duration_s + 60:
            continue

        rows.append(
            {
                "onset_s": onset_s,
                "duration_s": float(event_duration),
                "type": event_type,
                "apnea": event_type.upper().startswith(("HYP", "APNEA")),
            }
        )

    return pd.DataFrame(rows, columns=["onset_s", "duration_s", "type", "apnea"])


def _mark_event_epochs(events: pd.DataFrame, n_epochs: int, epoch_s: float) -> np.ndarray:
    labels = np.zeros(n_epochs, dtype=np.int8)
    if events is None or events.empty:
        return labels
    for row in events.itertuples(index=False):
        if not row.apnea:
            continue
        start = max(0, int(row.onset_s // epoch_s))
        end = min(n_epochs - 1, int((row.onset_s + row.duration_s - 1e-9) // epoch_s))
        labels[start : end + 1] = 1
    return labels


def _read_subject_details(root: Path) -> pd.DataFrame:
    path = _data_dir("ucddb", root) / "SubjectDetails.xls"
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_excel(path)
    df.columns = [str(c).strip() for c in df.columns]
    if "Study Number" in df.columns:
        df["Study Number"] = df["Study Number"].astype(str).str.strip()
    return df


def inspect_apnea_ecg(record_id: str, root: Path | None = None) -> RecordInfo:
    data_dir = _data_dir("apnea-ecg", root)
    base_id = record_id[:-2] if record_id.endswith("er") else record_id
    header = parse_wfdb_header(data_dir / f"{record_id}.hea")

    annotation_path = data_dir / record_id
    if not (data_dir / f"{record_id}.apn").exists():
        annotation_path = data_dir / base_id
    ann = wfdb.rdann(str(annotation_path), "apn")
    symbols = np.asarray(ann.symbol)
    labels = (symbols == "A").astype(np.int8)

    warnings = _check_spacing(ann.sample, ann.fs, 60.0, f"apnea-ecg/{record_id}")
    notes = "minute-level labels; respiration channels require base + r records"
    if warnings:
        notes += "; " + "; ".join(warnings)

    resp_names = [n for n in header["signal_names"] if n.startswith("Resp")]
    return RecordInfo(
        dataset="apnea-ecg",
        record_id=record_id,
        subject_id=base_id,
        n_signals=header["n_signals"],
        signal_names=header["signal_names"],
        ecg_channel="ECG" if "ECG" in header["signal_names"] else None,
        resp_channel=";".join(resp_names) if resp_names else None,
        fs_ecg=header["fs"],
        fs_resp=header["fs"] if resp_names else None,
        duration_s=header["duration_s"],
        label_source=".apn",
        label_unit_s=60.0,
        n_label_units=len(labels),
        n_apnea_units=int(labels.sum()),
        n_events=None,
        n_apnea_events=None,
        event_types={"apnea_minutes": int(labels.sum()), "normal_minutes": int((labels == 0).sum())},
        sleep_efficiency=None,
        sleep_fraction=None,
        notes=notes,
    )


def inspect_slpdb(
    record_id: str,
    root: Path | None = None,
    subject_details: pd.DataFrame | None = None,
) -> RecordInfo:
    data_dir = _data_dir("slpdb", root)
    header = parse_wfdb_header(data_dir / f"{record_id}.hea")
    ann = wfdb.rdann(str(data_dir / record_id), "st")

    stages = []
    event_types: dict[str, int] = {}
    apnea_units = 0
    n_events = 0
    n_apnea_events = 0

    for aux in ann.aux_note:
        tokens = _clean_aux(aux).split()
        if not tokens:
            stages.append("")
            continue
        stages.append(tokens[0])
        events = tokens[1:]
        n_events += len(events)
        has_apnea = False
        for token in events:
            event_types[token] = event_types.get(token, 0) + 1
            if token in APNEA_EVENTS:
                n_apnea_events += 1
                has_apnea = True
        if has_apnea:
            apnea_units += 1

    stages_arr = np.asarray(stages, dtype=object)
    sleep_units = int(np.isin(stages_arr, list(SLEEP_STAGES)).sum())
    sleep_fraction = sleep_units / len(stages_arr) if len(stages_arr) else None

    resp_names = [n for n in header["signal_names"] if n.startswith("Resp")]
    resp_channel = next((n for n in SLPDB_RESP_PREFERENCE if n in resp_names), None)
    if resp_channel is None and resp_names:
        resp_channel = resp_names[0]

    subject_match = re.match(r"(slp\d+)", record_id)
    subject_id = subject_match.group(1) if subject_match else record_id

    warnings = _check_spacing(ann.sample, ann.fs, 30.0, f"slpdb/{record_id}")
    notes = "30-s epoch labels (stage + event tokens); no event durations"
    if warnings:
        notes += "; " + "; ".join(warnings)

    return RecordInfo(
        dataset="slpdb",
        record_id=record_id,
        subject_id=subject_id,
        n_signals=header["n_signals"],
        signal_names=header["signal_names"],
        ecg_channel="ECG" if "ECG" in header["signal_names"] else None,
        resp_channel=resp_channel,
        fs_ecg=header["fs"],
        fs_resp=header["fs"] if resp_names else None,
        duration_s=header["duration_s"],
        label_source=".st",
        label_unit_s=30.0,
        n_label_units=len(stages_arr),
        n_apnea_units=apnea_units,
        n_events=n_events,
        n_apnea_events=n_apnea_events,
        event_types=event_types,
        sleep_efficiency=None,
        sleep_fraction=sleep_fraction,
        notes=notes,
    )


def inspect_ucddb(
    record_id: str,
    root: Path | None = None,
    subject_details: pd.DataFrame | None = None,
) -> RecordInfo:
    data_dir = _data_dir("ucddb", root)
    rec_path = data_dir / f"{record_id}.rec"
    stage_path = data_dir / f"{record_id}_stage.txt"
    event_path = data_dir / f"{record_id}_respevt.txt"

    reader = pyedflib.EdfReader(str(rec_path))
    try:
        signal_names = list(reader.getSignalLabels())
        fs_by_name = {
            name: float(reader.getSampleFrequency(i))
            for i, name in enumerate(signal_names)
        }
        duration_s = float(reader.getFileDuration())
        start_datetime = reader.getStartdatetime()
    finally:
        reader.close()

    stage = _read_ucddb_stage(stage_path)
    events = parse_ucddb_events(event_path, start_datetime, duration_s)
    epoch_s = 30.0
    labels = _mark_event_epochs(events, len(stage), epoch_s)

    sleep_fraction = (
        float(np.isin(stage, [1, 2, 3, 4, 5]).mean()) if len(stage) else None
    )
    event_types: dict[str, int] = {}
    for event_type in events["type"]:
        event_types[event_type] = event_types.get(event_type, 0) + 1

    resp_channel = next((n for n in UCDDB_RESP_PREFERENCE if n in signal_names), None)
    ahi = None
    sleep_efficiency = None
    if subject_details is not None and not subject_details.empty:
        match = subject_details[
            subject_details["Study Number"].str.upper() == record_id.upper()
        ]
        if not match.empty:
            if "PSG AHI" in match.columns:
                ahi = float(match.iloc[0]["PSG AHI"])
            if "Sleep Efficiency (%)" in match.columns:
                sleep_efficiency = float(match.iloc[0]["Sleep Efficiency (%)"])

    warnings = []
    expected_epochs = duration_s / epoch_s
    if abs(len(stage) - expected_epochs) > 2:
        warnings.append(
            f"stage epochs {len(stage)} != duration/30 {expected_epochs:.1f}"
        )
    notes = "use .rec (PSG), not _lifecard.edf (ECG-only); events have onset+duration"
    if warnings:
        notes += "; " + "; ".join(warnings)

    return RecordInfo(
        dataset="ucddb",
        record_id=record_id,
        subject_id=record_id,
        n_signals=len(signal_names),
        signal_names=signal_names,
        ecg_channel="ECG" if "ECG" in signal_names else None,
        resp_channel=resp_channel,
        fs_ecg=fs_by_name.get("ECG"),
        fs_resp=fs_by_name.get(resp_channel) if resp_channel else None,
        duration_s=duration_s,
        label_source="_stage.txt + _respevt.txt",
        label_unit_s=epoch_s,
        n_label_units=len(stage),
        n_apnea_units=int(labels.sum()),
        n_events=len(events),
        n_apnea_events=int(events["apnea"].sum()) if not events.empty else 0,
        event_types=event_types,
        ahi=ahi,
        sleep_efficiency=sleep_efficiency,
        sleep_fraction=sleep_fraction,
        notes=notes,
    )


def inspect_record(
    dataset: str,
    record_id: str,
    root: Path | None = None,
    subject_details: pd.DataFrame | None = None,
) -> RecordInfo:
    if dataset == "apnea-ecg":
        return inspect_apnea_ecg(record_id, root)
    if dataset == "slpdb":
        return inspect_slpdb(record_id, root, subject_details)
    if dataset == "ucddb":
        if subject_details is None:
            subject_details = _read_subject_details(Path(root) if root else DEFAULT_ROOT)
        return inspect_ucddb(record_id, root, subject_details)
    raise ValueError(f"Unknown dataset: {dataset}")


def list_records(dataset: str, root: Path | None = None) -> list[str]:
    data_dir = _data_dir(dataset, root)
    if dataset == "apnea-ecg":
        return sorted(path.stem for path in data_dir.glob("*er.hea"))
    if dataset == "slpdb":
        records = [
            line.strip()
            for line in (data_dir / "RECORDS").read_text().splitlines()
            if line.strip()
        ]
        return [r for r in records if (data_dir / f"{r}.hea").exists()]
    if dataset == "ucddb":
        records = []
        for line in (data_dir / "RECORDS").read_text().splitlines():
            record = line.strip()
            if not record or record.endswith(".edf"):
                continue
            if record.endswith(".rec"):
                record = record[:-4]
            records.append(record)
        return [r for r in records if (data_dir / f"{r}.rec").exists()]
    raise ValueError(f"Unknown dataset: {dataset}")


def _load_apnea_ecg(record_id: str, info: RecordInfo, root: Path | None = None) -> RecordData:
    data_dir = _data_dir("apnea-ecg", root)
    base_id = record_id[:-2] if record_id.endswith("er") else record_id

    ecg_record = wfdb.rdrecord(str(data_dir / base_id))
    resp_record = wfdb.rdrecord(str(data_dir / f"{base_id}r"))

    ecg = ecg_record.p_signal[:, 0].astype(np.float32)
    resp_names = list(resp_record.sig_name)
    resp_indices = [i for i, name in enumerate(resp_names) if name.startswith("Resp")]
    resp = resp_record.p_signal[:, resp_indices].astype(np.float32)

    n = min(len(ecg), len(resp))
    ecg = ecg[:n]
    resp = resp[:n]

    ann = wfdb.rdann(str(data_dir / record_id), "apn")
    labels = (np.asarray(ann.symbol) == "A").astype(np.int8)

    return RecordData(
        info=info,
        signals={"ecg": ecg, "resp": resp},
        sampling_rates={"ecg": float(ecg_record.fs), "resp": float(resp_record.fs)},
        label_times_s=ann.sample / ann.fs,
        labels=labels,
        sleep_stage=None,
        events=None,
    )


def _load_slpdb(record_id: str, info: RecordInfo, root: Path | None = None) -> RecordData:
    data_dir = _data_dir("slpdb", root)
    record = wfdb.rdrecord(str(data_dir / record_id))
    names = list(record.sig_name)

    ecg_idx = names.index("ECG")
    resp_idx = names.index(info.resp_channel) if info.resp_channel else None
    if resp_idx is None:
        raise ValueError(f"No respiratory channel found in slpdb/{record_id}")

    ecg = record.p_signal[:, ecg_idx].astype(np.float32)
    resp = record.p_signal[:, resp_idx].astype(np.float32)

    ann = wfdb.rdann(str(data_dir / record_id), "st")
    stages = []
    labels = []
    event_rows = []
    for i, aux in enumerate(ann.aux_note):
        tokens = _clean_aux(aux).split()
        stage = tokens[0] if tokens else ""
        events = tokens[1:]
        stages.append(stage)
        labels.append(int(any(t in APNEA_EVENTS for t in events)))
        event_rows.append(
            {
                "epoch": i,
                "onset_s": ann.sample[i] / ann.fs,
                "stage": stage,
                "events": ";".join(events),
            }
        )

    return RecordData(
        info=info,
        signals={"ecg": ecg, "resp": resp},
        sampling_rates={"ecg": float(record.fs), "resp": float(record.fs)},
        label_times_s=ann.sample / ann.fs,
        labels=np.asarray(labels, dtype=np.int8),
        sleep_stage=np.asarray(stages, dtype=object),
        events=pd.DataFrame(event_rows),
    )


def _load_ucddb(record_id: str, info: RecordInfo, root: Path | None = None) -> RecordData:
    data_dir = _data_dir("ucddb", root)
    rec_path = data_dir / f"{record_id}.rec"
    stage_path = data_dir / f"{record_id}_stage.txt"
    event_path = data_dir / f"{record_id}_respevt.txt"

    reader = pyedflib.EdfReader(str(rec_path))
    try:
        signal_names = list(reader.getSignalLabels())
        ecg_idx = signal_names.index("ECG")
        resp_idx = signal_names.index(info.resp_channel) if info.resp_channel else None
        if resp_idx is None:
            raise ValueError(f"No respiratory channel found in ucddb/{record_id}")

        ecg = reader.readSignal(ecg_idx).astype(np.float32)
        resp = reader.readSignal(resp_idx).astype(np.float32)
        fs_ecg = float(reader.getSampleFrequency(ecg_idx))
        fs_resp = float(reader.getSampleFrequency(resp_idx))
        start_datetime = reader.getStartdatetime()
        duration_s = float(reader.getFileDuration())
    finally:
        reader.close()

    stage = _read_ucddb_stage(stage_path)
    events = parse_ucddb_events(event_path, start_datetime, duration_s)
    labels = _mark_event_epochs(events, len(stage), 30.0)

    return RecordData(
        info=info,
        signals={"ecg": ecg, "resp": resp},
        sampling_rates={"ecg": fs_ecg, "resp": fs_resp},
        label_times_s=np.arange(len(stage), dtype=float) * 30.0,
        labels=labels,
        sleep_stage=stage,
        events=events,
    )


def load_record(dataset: str, record_id: str, root: Path | None = None) -> RecordData:
    subject_details = None
    if dataset == "ucddb":
        subject_details = _read_subject_details(Path(root) if root else DEFAULT_ROOT)
    info = inspect_record(dataset, record_id, root, subject_details)
    if dataset == "apnea-ecg":
        return _load_apnea_ecg(record_id, info, root)
    if dataset == "slpdb":
        return _load_slpdb(record_id, info, root)
    if dataset == "ucddb":
        return _load_ucddb(record_id, info, root)
    raise ValueError(f"Unknown dataset: {dataset}")


def load_labels(dataset: str, record_id: str, root: Path | None = None) -> RecordData:
    subject_details = None
    if dataset == "ucddb":
        subject_details = _read_subject_details(Path(root) if root else DEFAULT_ROOT)
    info = inspect_record(dataset, record_id, root, subject_details)

    if dataset == "apnea-ecg":
        data_dir = _data_dir(dataset, root)
        ann = wfdb.rdann(str(data_dir / record_id), "apn")
        labels = (np.asarray(ann.symbol) == "A").astype(np.int8)
        return RecordData(
            info=info,
            signals={},
            sampling_rates={},
            label_times_s=ann.sample / ann.fs,
            labels=labels,
        )

    if dataset == "slpdb":
        data_dir = _data_dir(dataset, root)
        ann = wfdb.rdann(str(data_dir / record_id), "st")
        stages = []
        labels = []
        event_rows = []
        for i, aux in enumerate(ann.aux_note):
            tokens = _clean_aux(aux).split()
            stage = tokens[0] if tokens else ""
            events = tokens[1:]
            stages.append(stage)
            labels.append(int(any(t in APNEA_EVENTS for t in events)))
            event_rows.append(
                {
                    "epoch": i,
                    "onset_s": ann.sample[i] / ann.fs,
                    "stage": stage,
                    "events": ";".join(events),
                }
            )
        return RecordData(
            info=info,
            signals={},
            sampling_rates={},
            label_times_s=ann.sample / ann.fs,
            labels=np.asarray(labels, dtype=np.int8),
            sleep_stage=np.asarray(stages, dtype=object),
            events=pd.DataFrame(event_rows),
        )

    if dataset == "ucddb":
        data_dir = _data_dir(dataset, root)
        reader = pyedflib.EdfReader(str(data_dir / f"{record_id}.rec"))
        try:
            start_datetime = reader.getStartdatetime()
            duration_s = float(reader.getFileDuration())
        finally:
            reader.close()
        stage = _read_ucddb_stage(data_dir / f"{record_id}_stage.txt")
        events = parse_ucddb_events(
            data_dir / f"{record_id}_respevt.txt", start_datetime, duration_s
        )
        labels = _mark_event_epochs(events, len(stage), 30.0)
        return RecordData(
            info=info,
            signals={},
            sampling_rates={},
            label_times_s=np.arange(len(stage), dtype=float) * 30.0,
            labels=labels,
            sleep_stage=stage,
            events=events,
        )

    raise ValueError(f"Unknown dataset: {dataset}")


def info_to_row(info: RecordInfo) -> dict[str, Any]:
    return {
        "dataset": info.dataset,
        "record_id": info.record_id,
        "subject_id": info.subject_id,
        "n_signals": info.n_signals,
        "signal_names": ";".join(info.signal_names),
        "ecg_channel": info.ecg_channel,
        "resp_channel": info.resp_channel,
        "fs_ecg": info.fs_ecg,
        "fs_resp": info.fs_resp,
        "duration_s": round(info.duration_s, 1),
        "duration_h": round(info.duration_s / 3600, 2),
        "label_source": info.label_source,
        "label_unit_s": info.label_unit_s,
        "n_label_units": info.n_label_units,
        "n_apnea_units": info.n_apnea_units,
        "apnea_prevalence": round(info.apnea_prevalence, 4),
        "n_events": info.n_events,
        "n_apnea_events": info.n_apnea_events,
        "event_types": json.dumps(info.event_types, sort_keys=True),
        "ahi": info.ahi,
        "sleep_efficiency": info.sleep_efficiency,
        "sleep_fraction": round(info.sleep_fraction, 4)
        if info.sleep_fraction is not None
        else None,
        "notes": info.notes,
    }
