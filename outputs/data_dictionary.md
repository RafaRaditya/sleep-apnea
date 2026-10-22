# Data Dictionary — `dataset_inventory.csv` & `dataset_summary.json`

## Hubungan dua file

| File | Isi | Fungsi |
|---|---|---|
| `dataset_inventory.csv` | 51 baris, satu baris per rekaman | Buku besar / bukti audit per record |
| `dataset_summary.json` | Agregat per dataset + overall | Angka utama untuk paper/abstrak |

Keduanya dihasilkan oleh `code/dataset_inventory.py` langsung dari file mentah WFDB/EDF. CSV = detail; JSON = ringkasan.

## Kolom `dataset_inventory.csv`

### Identitas
| Kolom | Arti |
|---|---|
| `dataset` | `apnea-ecg`, `slpdb`, atau `ucddb` |
| `record_id` | Nama record, mis. `a01er`, `slp01a`, `ucddb002` |
| `subject_id` | ID subject; `slp01a` dan `slp01b` sama-sama `slp01` |

### Sinyal
| Kolom | Arti |
|---|---|
| `n_signals` | Jumlah kanal dalam record |
| `signal_names` | Daftar kanal, dipisah `;` |
| `ecg_channel` | Nama kanal ECG yang dipakai |
| `resp_channel` | Kanal respirasi yang dipakai |
| `fs_ecg`, `fs_resp` | Sampling rate ECG dan respirasi (Hz) |
| `duration_s`, `duration_h` | Durasi rekaman |

### Label
| Kolom | Arti |
|---|---|
| `label_source` | File label: `.apn`, `.st`, atau `_stage.txt + _respevt.txt` |
| `label_unit_s` | Panjang satu unit label: 60 detik, 30 detik |
| `n_label_units` | Jumlah unit label (menit/epoch) |
| `n_apnea_units` | Jumlah unit yang ditandai apnea |
| `apnea_prevalence` | `n_apnea_units / n_label_units` |

### Event & klinis
| Kolom | Arti |
|---|---|
| `n_events` | Semua event yang dianotasi; kosong untuk Apnea-ECG |
| `n_apnea_events` | Event apnea/hypopnea saja |
| `event_types` | Rincian tipe event dalam format JSON |
| `ahi` | AHI klinis; hanya UCDDB (dari `SubjectDetails.xls`) |
| `sleep_efficiency` | Sleep efficiency klinis; hanya UCDDB |
| `sleep_fraction` | Fraksi epoch yang dinilai tidur (bukan sleep efficiency) |
| `notes` | Catatan teknis per record |

## Field `dataset_summary.json`

| Field | Arti |
|---|---|
| `n_records` | Jumlah record |
| `n_subjects` | Jumlah subject unik |
| `total_duration_h` | Total durasi rekaman (jam) |
| `n_label_units` | Total unit label |
| `n_apnea_units` | Total unit apnea |
| `apnea_prevalence` | Total unit apnea / total unit label |
| `n_events` | Total event (null untuk Apnea-ECG) |
| `n_apnea_events` | Total event apnea/hypopnea |
| `mean_ahi` | Rata-rata AHI (hanya UCDDB) |
| `mean_sleep_fraction` | Rata-rata fraksi epoch tidur |
| `warnings` | Daftar masalah validasi; `[]` berarti semua record lolos cek dasar |

## Cara membaca tiga contoh baris

1. **`apnea-ecg / a01er`**
   - 5 kanal: ECG, Resp C, Resp A, Resp N, SpO2.
   - 489 unit label (menit), 470 apnea → prevalensi 96,1%.
   - `n_events` kosong karena label hanya per menit, tidak ada onset/durasi event.

2. **`slpdb / slp01a`**
   - 4 kanal: ECG, BP, EEG, Resp (sum).
   - 240 epoch (30 detik), 37 apnea → 15,4%.
   - 38 event apnea; `event_types` juga mencatat `L/LA` (leg movement) yang bukan apnea.

3. **`ucddb / ucddb002`**
   - 14 kanal PSG: ECG 128 Hz; Flow/ribcage/abdo 8 Hz.
   - 748 epoch (30 detik), 182 apnea → 24,3%.
   - 124 event, AHI 23, sleep efficiency 84%.

## Pertanyaan yang sering muncul

**Kenapa 51 record tapi 49 subject?**
Karena `slp01a`/`slp01b` dan `slp02a`/`slp02b` adalah subject yang sama. Jumlah subject: 8 + 16 + 25 = 49.

**Kenapa Apnea-ECG cuma 8?**
Hanya record `*er` yang punya kanal respirasi (`Resp C/A/N`). Record lain hanya ECG.

**Kenapa `n_events` Apnea-ECG kosong?**
Label `.apn` hanya memberi A/N per menit, tanpa onset/durasi event.

**Kenapa prevalensi UCDDB 21,2% tapi AHI 24,2?**
Beda satuan. Prevalensi = fraksi epoch 30 detik yang beririsan dengan event. AHI = jumlah event per jam tidur dari skoring klinis. Keduanya tidak bisa dibandingkan langsung.

**Apa beda `n_events` dan `n_apnea_events`?**
`n_events` menghitung semua anotasi (termasuk `POSSIBLE`, leg movement). `n_apnea_events` hanya `HYP-*`/`APNEA-*` atau `H/HA/OA/X/CA/CAA`.

**Apa beda `sleep_fraction` dan `sleep_efficiency`?**
`sleep_fraction` = proporsi epoch yang dinilai tidur. `sleep_efficiency` = metrik klinis UCDDB (waktu tidur / waktu di tempat tidur). Jangan disamakan.

**Apa arti `warnings: []`?**
Semua 51 record lolos pengecekan dasar: ada ECG, ada respirasi, ada label, durasi positif. Kalau ada masalah, warning akan muncul di sini.

## Cara presentasi 1 menit

1. Buka JSON, sebut 4 angka: **51 record, 49 subject, 325,08 jam, 26,79% prevalensi total**.
2. Buka CSV, filter satu dataset, tunjukkan satu baris contoh (mis. `a01er` dan `ucddb002`).
3. Jelaskan kolom kunci: `resp_channel`, `label_unit_s`, `n_apnea_units`, `apnea_prevalence`.
4. Tutup dengan caveat: **CSV adalah bukti per record; JSON adalah ringkasan. Keduanya dihasilkan dari file mentah, bukan dari dokumentasi dataset.**
