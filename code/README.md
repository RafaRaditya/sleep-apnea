# Sleep Apnea Research — Inventory, Harmonization & Loaders

Script ini memvalidasi 51 rekaman dari tiga dataset (Apnea-ECG, MIT-BIH PSG/slpdb, UCDDB), menyediakan loader standar ECG + respirasi + label, dan menghasilkan dataset harmonized pada 100 Hz.

## Cara pakai

```bash
pip install -r requirements.txt

# Bangun inventory 51 record -> outputs/dataset_inventory.csv + dataset_summary.json
python code/dataset_inventory.py

# Uji loader pada satu record per dataset
python code/dataset_inventory.py --self-test

# EDA: window 30/10 detik, sleep vs wake, durasi event, plot contoh
python code/eda.py

# Verifikasi independen: hitung ulang dari file mentah + checksum resmi PhysioNet
python code/verify_inventory.py --checksums all

# Harmonisasi: 100 Hz, ECG + respirasi, label 30 detik -> outputs/harmonized/
python code/harmonize.py

# Verifikasi hasil harmonisasi
python code/verify_harmonized.py
```

Output verifikasi inventory:

- `outputs/verification_report.md` — ringkasan PASS/FAIL.
- `outputs/verification_results.csv` — status per record.

Output utama:

- `outputs/dataset_inventory.csv` — satu baris per record: sinyal, fs, durasi, sumber label, jumlah unit apnea, prevalensi, event, AHI.
- `outputs/dataset_summary.json` — ringkasan per dataset dan keseluruhan.
- `outputs/eda_report.md` — hasil EDA dan implikasinya.
- `outputs/figures/` — plot prevalensi per record, sleep fraction, durasi event, dan contoh window apnea vs normal.
- `outputs/harmonized/` — sinyal 100 Hz + label 30 detik + manifest; detail format di `outputs/harmonized/README.md`.

## Format harmonized

Setiap record disimpan sebagai `outputs/harmonized/signals/{dataset}__{record_id}.npz`:

| Key | Arti |
|---|---|
| `ecg`, `resp` | Sinyal float32 pada 100 Hz |
| `labels_30s` | Label apnea 0/1 per epoch 30 detik |
| `label_mask_30s` | 1 = epoch berlabel; 0 = unlabeled (drop saat training) |
| `sleep_stage_30s` | Kode stage: -1 unknown, 0 wake, 1 N1, 2 N2, 3 N3, 4 N4, 5 REM, 6 movement |
| `native_labels`, `native_label_times_s`, `native_label_unit_s` | Label asli (60/30 detik) |
| `event_onset_s`, `event_duration_s`, `event_apnea`, `event_type` | Event-level UCDDB |
| `slpdb_event_epoch`, `slpdb_event_tokens` | Event per epoch slpdb |
| `meta_json` | Metadata lengkap |

Sinyal di harmonized **belum difilter dan belum dinormalisasi**. Filtering dan normalisasi/clipping dilakukan setelah split per subject (Tahap 5) agar tidak ada leakage.

## API loader

```python
from code.loaders import load_record, list_records

data = load_record("ucddb", "ucddb002")
data.signals["ecg"]        # 1D ECG
data.signals["resp"]       # 1D respiratory channel
data.sampling_rates        # {"ecg": 128.0, "resp": 8.0}
data.labels                # 0/1 per label unit
data.label_times_s         # waktu mulai tiap label unit
data.sleep_stage           # None untuk Apnea-ECG
data.events                # DataFrame onset/duration untuk UCDDB
```

## Keputusan dataset yang sudah ditangani

- **Apnea-ECG**: hanya 8 record `*er` yang punya respirasi (`Resp C/A/N`) + SpO2. Header `a01er.hea` adalah multi-file; loader membaca `a01` (ECG) dan `a01r` (respirasi) terpisah lalu memotong ke panjang minimum.
- **slpdb**: label `.st` adalah epoch 30 detik; token pertama = sleep stage, sisanya event. Event apnea = `H/HA/OA/X/CA/CAA`; `L/LA/A/MT` bukan apnea. `slp41` dan `slp45` tidak punya event apnea.
- **UCDDB**: memakai `.rec` (PSG EDF 14 kanal), bukan `_lifecard.edf` (3 kanal ECG saja). `_stage.txt` adalah epoch 30 detik; `_respevt.txt` berisi onset + durasi event. Waktu event dikonversi memakai start datetime EDF dan tahan lintas tengah malam.

## Batasan yang perlu diingat

- Label Apnea-ECG per **menit**, slpdb per **30 detik**, UCDDB punya onset + durasi. Ini memengaruhi validitas eksperimen window 10 detik.
- `slp01a`/`slp01b` dan `slp02a`/`slp02b` adalah subject yang sama; split harus per subject, bukan per record.
- UCDDB `Flow` 8 Hz; ECG 128 Hz. Harmonisasi ke 100 Hz sudah dilakukan dan disimpan di `outputs/harmonized/`.
