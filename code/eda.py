from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from loaders import (
    DEFAULT_ROOT,
    list_records,
    load_labels,
    load_record,
)


DATASETS = ["apnea-ecg", "slpdb", "ucddb"]
COLORS = {"apnea-ecg": "#d62728", "slpdb": "#1f77b4", "ucddb": "#2ca02c"}
SLEEP_STAGES_SLPDB = {"1", "2", "3", "4", "R"}
SLEEP_STAGES_UCDDB = {1, 2, 3, 4, 5}


def mark_event_epochs(events: pd.DataFrame, n_epochs: int, epoch_s: float) -> np.ndarray:
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


def build_window_table(root: Path) -> pd.DataFrame:
    rows = []
    for dataset in DATASETS:
        for record_id in list_records(dataset, root):
            data = load_labels(dataset, record_id, root)
            labels = data.labels

            agreement = None
            propagated_fp = None
            propagated_fn = None

            if dataset == "apnea-ecg":
                labels_30 = np.repeat(labels, 2)
                labels_10 = np.repeat(labels, 6)
            elif dataset == "slpdb":
                labels_30 = labels
                labels_10 = np.repeat(labels, 3)
            else:
                labels_30 = labels
                labels_10 = mark_event_epochs(data.events, len(labels) * 3, 10.0)
                parent_10 = np.repeat(labels_30, 3)
                agreement = float((parent_10 == labels_10).mean())
                propagated_fp = int(((parent_10 == 1) & (labels_10 == 0)).sum())
                propagated_fn = int(((parent_10 == 0) & (labels_10 == 1)).sum())

            rows.append(
                {
                    "dataset": dataset,
                    "record_id": record_id,
                    "subject_id": data.info.subject_id,
                    "n30": len(labels_30),
                    "pos30": int(labels_30.sum()),
                    "prev30": float(labels_30.mean()) if len(labels_30) else 0.0,
                    "n10": len(labels_10),
                    "pos10": int(labels_10.sum()),
                    "prev10": float(labels_10.mean()) if len(labels_10) else 0.0,
                    "propagated_agreement_10": agreement,
                    "propagated_fp_10": propagated_fp,
                    "propagated_fn_10": propagated_fn,
                }
            )
    return pd.DataFrame(rows)


def summarize_windows(window_df: pd.DataFrame) -> pd.DataFrame:
    agg = (
        window_df.groupby("dataset", sort=False)
        .agg(
            n_records=("record_id", "count"),
            n30=("n30", "sum"),
            pos30=("pos30", "sum"),
            n10=("n10", "sum"),
            pos10=("pos10", "sum"),
        )
        .reset_index()
    )
    agg["prev30"] = (agg["pos30"] / agg["n30"]).round(4)
    agg["prev10"] = (agg["pos10"] / agg["n10"]).round(4)
    return agg


def sleep_wake_stats(root: Path) -> pd.DataFrame:
    rows = []
    for dataset in ["slpdb", "ucddb"]:
        for record_id in list_records(dataset, root):
            data = load_labels(dataset, record_id, root)
            stage = data.sleep_stage
            if stage is None:
                continue
            if dataset == "slpdb":
                sleep_mask = np.isin(stage, list(SLEEP_STAGES_SLPDB))
            else:
                sleep_mask = np.isin(stage, list(SLEEP_STAGES_UCDDB))
            labels = data.labels

            def rate(mask: np.ndarray) -> float | None:
                return float(labels[mask].mean()) if mask.any() else None

            rows.append(
                {
                    "dataset": dataset,
                    "record_id": record_id,
                    "n_epochs": len(stage),
                    "sleep_fraction": float(sleep_mask.mean()),
                    "apnea_prev_sleep": rate(sleep_mask),
                    "apnea_prev_wake": rate(~sleep_mask),
                }
            )
    return pd.DataFrame(rows)


def event_duration_stats(root: Path) -> pd.DataFrame:
    rows = []
    for record_id in list_records("ucddb", root):
        data = load_labels("ucddb", record_id, root)
        if data.events is None or data.events.empty:
            continue
        events = data.events[data.events["apnea"]]
        for row in events.itertuples(index=False):
            rows.append(
                {
                    "record_id": record_id,
                    "duration_s": row.duration_s,
                    "type": row.type,
                }
            )
    return pd.DataFrame(rows)


def plot_prevalence_by_record(inventory: pd.DataFrame, out_path: Path) -> None:
    df = inventory.sort_values(["dataset", "apnea_prevalence"], ascending=[True, False])
    colors = [COLORS[d] for d in df["dataset"]]

    fig, ax = plt.subplots(figsize=(10, max(7, len(df) * 0.24)))
    ax.barh(range(len(df)), df["apnea_prevalence"], color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["record_id"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Apnea prevalence (label units)")
    ax.set_title("Label prevalence per recording")
    overall = inventory["n_apnea_units"].sum() / inventory["n_label_units"].sum()
    ax.axvline(overall, color="black", linestyle="--", linewidth=1, label=f"overall={overall:.2f}")
    handles = [
        plt.Line2D([0], [0], color=color, lw=6, label=dataset)
        for dataset, color in COLORS.items()
    ]
    ax.legend(handles=handles + [plt.Line2D([0], [0], color="black", lw=1, ls="--", label="overall")], fontsize=7)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_sleep_fraction(inventory: pd.DataFrame, out_path: Path) -> None:
    df = inventory.dropna(subset=["sleep_fraction"]).sort_values(
        ["dataset", "sleep_fraction"], ascending=[True, False]
    )
    colors = [COLORS[d] for d in df["dataset"]]
    fig, ax = plt.subplots(figsize=(10, max(5, len(df) * 0.24)))
    ax.barh(range(len(df)), df["sleep_fraction"], color=colors)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["record_id"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("Sleep epochs / all epochs")
    ax.set_title("Sleep fraction per recording")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_event_durations(durations: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(durations["duration_s"], bins=np.arange(0, 70, 2), color="#2ca02c", edgecolor="black")
    median = float(durations["duration_s"].median())
    ax.axvline(median, color="black", linestyle="--", label=f"median={median:.1f}s")
    ax.set_xlabel("Event duration (s)")
    ax.set_ylabel("Count")
    ax.set_title("UCDDB apnea/hypopnea event durations")
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def plot_example_windows(root: Path, out_path: Path) -> None:
    examples = {"apnea-ecg": "a03er", "slpdb": "slp03", "ucddb": "ucddb002"}
    fig, axes = plt.subplots(len(examples), 2, figsize=(13, 9), sharex=False)

    for row, (dataset, record_id) in enumerate(examples.items()):
        data = load_record(dataset, record_id, root)
        labels = data.labels
        pos = np.where(labels == 1)[0]
        neg = np.where(labels == 0)[0]
        if len(pos) == 0 or len(neg) == 0:
            continue

        for col, (idx, condition) in enumerate([(pos[0], "apnea"), (neg[0], "normal")]):
            t0 = float(data.label_times_s[idx])
            window_s = 30.0
            ecg = data.signals["ecg"]
            resp = np.asarray(data.signals["resp"])
            fs_ecg = data.sampling_rates["ecg"]
            fs_resp = data.sampling_rates["resp"]

            i0 = int(t0 * fs_ecg)
            i1 = min(len(ecg), i0 + int(window_s * fs_ecg))
            j0 = int(t0 * fs_resp)
            j1 = min(len(resp), j0 + int(window_s * fs_resp))
            if i1 <= i0 or j1 <= j0:
                continue

            ax = axes[row, col]
            ax.plot(np.arange(i0, i1) / fs_ecg - t0, ecg[i0:i1], color="#1f77b4", lw=0.8)
            ax.set_ylabel("ECG", color="#1f77b4")
            ax.tick_params(axis="y", labelcolor="#1f77b4")

            resp_window = resp[j0:j1]
            if resp_window.ndim == 1:
                resp_window = resp_window[:, None]
            channel = resp_window.shape[1] - 1 if dataset == "apnea-ecg" else 0
            ax2 = ax.twinx()
            ax2.plot(
                np.arange(j0, j1) / fs_resp - t0,
                resp_window[:, channel],
                color="#ff7f0e",
                lw=0.8,
            )
            ax2.set_ylabel("Resp", color="#ff7f0e")
            ax2.tick_params(axis="y", labelcolor="#ff7f0e")

            ax.set_title(f"{dataset}/{record_id} | {condition} | t={t0:.0f}s")
            ax.set_xlabel("seconds")

    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def write_report(
    path: Path,
    inventory: pd.DataFrame,
    summary: dict,
    window_summary: pd.DataFrame,
    window_df: pd.DataFrame,
    sleep_df: pd.DataFrame,
    durations: pd.DataFrame,
) -> None:
    lines = []
    lines.append("# EDA Report — Sleep Apnea Multi-Dataset\n")
    lines.append("## Dataset summary\n")
    lines.append("```")
    lines.append(json.dumps({k: v for k, v in summary.items() if k != "warnings"}, indent=2))
    lines.append("```\n")

    lines.append("## Window estimates\n")
    lines.append("```")
    lines.append(window_summary.to_string(index=False))
    lines.append("```\n")

    ucddb = window_df[window_df["dataset"] == "ucddb"]
    if not ucddb.empty:
        total_10 = int(ucddb["n10"].sum())
        fp = int(ucddb["propagated_fp_10"].sum())
        fn = int(ucddb["propagated_fn_10"].sum())
        agreement = 1.0 - (fp + fn) / total_10 if total_10 else None
        lines.append("## UCDDB: 30-s label propagated to 10-s windows\n")
        lines.append(
            f"- 10-s windows: {total_10}\n"
            f"- Agreement with event-level 10-s labels: {agreement:.4f}\n"
            f"- Propagated positive but event negative: {fp}\n"
            f"- Propagated negative but event positive: {fn}\n"
        )

    lines.append("## Sleep vs wake apnea prevalence\n")
    if not sleep_df.empty:
        agg = (
            sleep_df.groupby("dataset", sort=False)
            .agg(
                n_records=("record_id", "count"),
                mean_sleep_fraction=("sleep_fraction", "mean"),
                mean_prev_sleep=("apnea_prev_sleep", "mean"),
                mean_prev_wake=("apnea_prev_wake", "mean"),
            )
            .reset_index()
        )
        lines.append("```")
        lines.append(agg.round(3).to_string(index=False))
        lines.append("```\n")

    lines.append("## UCDDB event duration\n")
    if not durations.empty:
        desc = durations["duration_s"].describe(percentiles=[0.1, 0.5, 0.9]).round(1)
        lines.append("```")
        lines.append(desc.to_string())
        lines.append("```")
        short = float((durations["duration_s"] < 10).mean())
        lines.append(f"\nFraction of events shorter than 10 s: {short:.3f}\n")

    lines.append("## Implications\n")
    lines.append("- **Window 10 s:** only UCDDB has event onsets and durations. For slpdb and Apnea-ECG, 10-s labels are inherited from 30-s/60-s units.")
    lines.append("- **Wake/sleep:** slpdb and UCDDB have sleep stages; Apnea-ECG does not. Excluding wake changes the label distribution and cannot be applied consistently to Apnea-ECG.")
    lines.append("- **Subject split:** slp01a/b and slp02a/b share subjects; group by subject, not record.")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="EDA for the sleep-apnea datasets.")
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()

    root = args.root.resolve()
    out_dir = args.out or root / "outputs"
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    inventory = pd.read_csv(out_dir / "dataset_inventory.csv")
    summary = json.loads((out_dir / "dataset_summary.json").read_text(encoding="utf-8"))

    print("Building window table...")
    window_df = build_window_table(root)
    window_summary = summarize_windows(window_df)
    window_df.to_csv(out_dir / "window_estimates_by_record.csv", index=False)

    print("Sleep/wake stats...")
    sleep_df = sleep_wake_stats(root)
    sleep_df.to_csv(out_dir / "sleep_wake_stats.csv", index=False)

    print("Event durations...")
    durations = event_duration_stats(root)
    durations.to_csv(out_dir / "ucddb_event_durations.csv", index=False)

    print("Plots...")
    plot_prevalence_by_record(inventory, fig_dir / "label_prevalence_by_record.png")
    plot_sleep_fraction(inventory, fig_dir / "sleep_fraction_by_record.png")
    if not durations.empty:
        plot_event_durations(durations, fig_dir / "ucddb_event_durations.png")
    plot_example_windows(root, fig_dir / "example_windows.png")

    write_report(
        out_dir / "eda_report.md",
        inventory,
        summary,
        window_summary,
        window_df,
        sleep_df,
        durations,
    )

    print("\nWindow estimates:")
    print(window_summary.to_string(index=False))
    print("\nWrote outputs/eda_report.md and outputs/figures/*.png")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
