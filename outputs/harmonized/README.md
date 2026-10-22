# Harmonized Dataset

Dihasilkan oleh `code/harmonize.py`. Tahap ini **hanya menyamakan format, sampling rate, kanal, dan label**. Sinyal masih **belum difilter dan belum dinormalisasi** — itu Tahap 5 (preprocessing) setelah split per subject, supaya tidak ada leakage.

## Isi folder

```
outputs/harmonized/
  manifest.csv                  # satu baris per record
  summary.json                  # agregat per dataset
  signals/
    {dataset}__{record_id}.npz  # sinyal + label + metadata array
    {dataset}__{record_id}.json # metadata (sidecar)
  verification_report.md        # hasil verify_harmonized.py
  verification_results.csv
```

## Cara membaca `manifest.csv`

Satu baris = satu record hasil harmonisasi. Kolom `native_*` menyimpan kondisi **asli** sebelum diproses; kolom `target_fs` dan turunannya adalah kondisi **setelah** harmonisasi.

| Kolom | Arti |
|---|---|
| `dataset` | `apnea-ecg`, `slpdb`, atau `ucddb` |
| `record_id` | Nama record (mis. `a01er`, `slp01a`, `ucddb002`) |
| `subject_id` | ID subject; `slp01a` dan `slp01b` sama-sama `slp01` |
| `native_fs_ecg` | Sampling rate ECG **asli** di file mentah (Hz) |
| `native_fs_resp` | Sampling rate respirasi **asli** di file mentah (Hz) |
| `target_fs` | Sampling rate hasil harmonisasi — **selalu 100.0** |
| `resp_channel` | Kanal respirasi yang dipilih (lihat aturan prioritas di bawah) |
| `duration_h` | Durasi sinyal setelah dipotong ke panjang minimum ECG/respirasi (jam) |
| `n_samples` | Jumlah sample per sinyal di `.npz`; selalu `duration_h × 3600 × target_fs` |
| `n_epochs_30s` | Jumlah epoch 30 detik yang utuh |
| `n_labeled_30s` | Epoch yang punya label (`label_mask_30s == 1`) |
| `n_apnea_30s` | Epoch berlabel apnea |
| `apnea_prevalence_30s` | `n_apnea_30s / n_labeled_30s` |
| `n_native_label_units` | Jumlah unit label asli (menit untuk Apnea-ECG, epoch 30 detik untuk slpdb/UCDDB) |
| `native_label_unit_s` | Panjang satu unit label asli: 60 atau 30 detik |
| `nan_fraction_ecg`, `nan_fraction_resp` | Fraksi NaN **sebelum** interpolasi |
| `max_nan_gap_resp_s` | Gap NaN terpanjang di respirasi **sebelum** interpolasi (detik) |
| `sleep_known_fraction` | Fraksi epoch yang punya sleep stage valid (bukan `-1`) |
| `n_sleep_30s` | Jumlah epoch yang dinilai tidur (stage 1–5); kosong untuk Apnea-ECG |
| `sleep_fraction_30s` | `n_sleep_30s / n_epochs_30s`; kosong untuk Apnea-ECG |
| `notes` | Catatan teknis dari loader |

### `native_fs_*` vs `target_fs`

`native_fs_ecg` dan `native_fs_resp` **sengaja berbeda-beda** — itu rate asli tiap dataset, disimpan untuk audit dan reproduksi:

```
slpdb,slp01a,...    native_fs_ecg=250.0  native_fs_resp=250.0  target_fs=100.0
ucddb,ucddb002,...  native_fs_ecg=128.0  native_fs_resp=8.0    target_fs=100.0
```

Sinyal yang tersimpan di `.npz` sudah **100 Hz** untuk semua record. Cek cepat:

```python
n_samples / (duration_h * 3600)  # harus sama dengan target_fs (100.0)
```

Contoh: `ucddb002` → `2247000 / (6.242 * 3600) = 100.0`.

### Kolom mana yang dipakai untuk apa

- **Cek kualitas / drop record**: `nan_fraction_*`, `max_nan_gap_resp_s`, `sleep_known_fraction`.
- **Hitung prevalensi**: `n_apnea_30s / n_labeled_30s` (pakai epoch berlabel saja, **bukan** `n_epochs_30s`).
- **Audit asal-usul**: `native_fs_*`, `n_native_label_units`, `native_label_unit_s`, `resp_channel`.
- **Validasi hasil harmonisasi**: `n_samples` harus konsisten dengan `duration_h × 3600 × target_fs`.

### Membaca dengan pandas

```python
import numpy as np
import pandas as pd

manifest = pd.read_csv("outputs/harmonized/manifest.csv")

# ringkasan per dataset
print(manifest.groupby("dataset")[["n_labeled_30s", "n_apnea_30s"]].sum())

# record dengan gap NaN respirasi terpanjang
print(manifest.nlargest(5, "max_nan_gap_resp_s")[["record_id", "max_nan_gap_resp_s"]])

# cek target rate konsisten
assert (manifest["target_fs"] == 100.0).all()
assert np.allclose(
    manifest["n_samples"] / (manifest["duration_h"] * 3600), 100.0, atol=0.05
)
```

## Isi `.npz`

| Key | Shape / dtype | Arti |
|---|---|---|
| `ecg` | `(n_samples,)` float32 | ECG pada 100 Hz |
| `resp` | `(n_samples,)` float32 | Satu kanal respirasi pada 100 Hz |
| `fs` | scalar float32 | 100.0 |
| `labels_30s` | `(n_epochs,)` int8 | Label apnea 0/1 per epoch 30 detik |
| `label_mask_30s` | `(n_epochs,)` int8 | 1 = epoch punya label; 0 = unlabeled |
| `label_times_30s` | `(n_epochs,)` float32 | Waktu mulai epoch (detik) |
| `sleep_stage_30s` | `(n_epochs,)` int8 | Kode sleep stage gabungan |
| `native_labels` | `(n_native,)` int8 | Label asli (menit/30 detik) |
| `native_label_times_s` | `(n_native,)` float32 | Waktu label asli |
| `native_label_unit_s` | scalar float32 | 60 (Apnea-ECG) atau 30 (slpdb/UCDDB) |
| `event_onset_s` | `(n_events,)` float32 | Khusus UCDDB |
| `event_duration_s` | `(n_events,)` float32 | Khusus UCDDB |
| `event_apnea` | `(n_events,)` int8 | Khusus UCDDB |
| `event_type` | `(n_events,)` unicode | Khusus UCDDB, mis. `HYP-O` |
| `slpdb_event_epoch` | `(n_native,)` int16 | Khusus slpdb |
| `slpdb_event_tokens` | `(n_native,)` unicode | Khusus slpdb, mis. `H A` |
| `meta_json` | scalar unicode | Metadata lengkap |

## Kode sleep stage

| Kode | Arti |
|---:|---|
| `-1` | Unknown / tidak tersedia |
| `0` | Wake |
| `1` | N1 |
| `2` | N2 |
| `3` | N3 |
| `4` | N4 |
| `5` | REM |
| `6` | Movement / other |

Mapping:
- slpdb: `W→0`, `1→1`, `2→2`, `3→3`, `4→4`, `R→5`, `MT/M→6`.
- UCDDB: `0→0`, `1→5` (REM), `2→1`, `3→2`, `4→3`, `5→4`, `6→6`, `7→-1`.
- Apnea-ECG: semua `-1` (tidak punya sleep stage).

## Aturan label

- **Apnea-ECG**: label `.apn` per 60 detik; setiap label menandai dua epoch 30 detik.
- **slpdb**: label `.st` per 30 detik; apnea jika ada token `H/HA/OA/X/CA/CAA`.
- **UCDDB**: label 30 detik dari irisan event `HYP-*`/`APNEA-*` dengan epoch.
- `label_mask_30s = 0` berarti epoch tidak punya label (mis. ujung rekaman yang tidak dianotasi). **Epoch ini harus di-drop saat training/evaluasi.**

## Pemilihan kanal

- ECG: kanal `ECG` di semua dataset.
- Respirasi:
  - Apnea-ECG: `Resp N` (nasal) dari record `*r`.
  - slpdb: prioritas `Resp (sum)` → `Resp (nasal)` → `Resp (chest)` → `Resp (abdominal)`; kanal aktual ada di `manifest.csv`.
  - UCDDB: `Flow`.
- SpO2 Apnea-ECG tidak dipakai pada tahap ini (di luar scope ECG + respirasi).

## Catatan penting

- `apnea_prevalence_30s` dihitung dari epoch 30 detik yang **punya label**. Karena Apnea-ECG dilipatgandakan dari menit ke 30 detik, prevalensi total (28,21%) berbeda dari prevalensi inventory berbasis unit asli (26,79%). Prevalensi per dataset tetap sama.
- `n_epochs_30s` bisa lebih besar dari `n_labeled_30s` jika anotasi berhenti sebelum sinyal berakhir. Gunakan `label_mask_30s`.
- Sinyal masih mentah (hanya resampling + interpolasi NaN). Filtering dan normalisasi/clipping dilakukan setelah split.

## Cara load

```python
import json
import numpy as np

path = "outputs/harmonized/signals/ucddb__ucddb002.npz"
data = np.load(path, allow_pickle=False)
ecg = data["ecg"]
resp = data["resp"]
labels = data["labels_30s"]
mask = data["label_mask_30s"]
sleep_stage = data["sleep_stage_30s"]
meta = json.loads(str(data["meta_json"].item()))
```
