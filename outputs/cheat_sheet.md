# Cheat Sheet — Cara Membaca `dataset_inventory.csv` & `dataset_summary.json`

## Model mental paling sederhana

- **CSV = tabel.** Satu baris = satu rekaman. Satu kolom = satu fakta tentang rekaman itu.
- **JSON = jumlah-jumlah dari tabel.** Tidak ada angka baru; semuanya agregasi dari CSV.

Kalau lupa, buka `outputs/data_dictionary.md`.

---

## 1. Membaca satu baris: `apnea-ecg / a01er`

| Kolom | Nilai | Artinya |
|---|---|---|
| `dataset` | `apnea-ecg` | Dari dataset Apnea-ECG |
| `record_id` | `a01er` | Nama rekaman |
| `subject_id` | `a01` | Subject-nya `a01` |
| `n_signals` | 5 | Ada 5 kanal |
| `signal_names` | `ECG;Resp C;Resp A;Resp N;SpO2` | ECG + 3 respirasi + SpO2 |
| `ecg_channel` | `ECG` | Kanal ECG yang dipakai |
| `resp_channel` | `Resp C;Resp A;Resp N` | Tiga kanal respirasi |
| `fs_ecg` | 100 | ECG 100 Hz |
| `fs_resp` | 100 | Respirasi 100 Hz |
| `duration_s` | 29570 | 29.570 detik |
| `duration_h` | 8.21 | 8,21 jam |
| `label_source` | `.apn` | Label dari file `.apn` |
| `label_unit_s` | 60 | Satu label = 1 menit |
| `n_label_units` | 489 | Ada 489 menit label |
| `n_apnea_units` | 470 | 470 menit ditandai apnea |
| `apnea_prevalence` | 0.9611 | 470/489 = 96,11% |
| `n_events` | *(kosong)* | Tidak ada onset/durasi event |
| `n_apnea_events` | *(kosong)* | Sama, tidak tersedia |
| `event_types` | `{"apnea_minutes":470,"normal_minutes":19}` | 470 apnea + 19 normal |
| `ahi` | *(kosong)* | Tidak dihitung dari file ini |
| `sleep_efficiency` | *(kosong)* | Apnea-ECG tidak punya sleep stage |
| `sleep_fraction` | *(kosong)* | Sama |
| `notes` | minute-level labels... | Catatan teknis |

**Kalimat satu baris:**
> “`a01er` adalah rekaman Apnea-ECG selama 8,21 jam dengan 5 kanal (ECG + tiga respirasi + SpO2), 489 label per menit, dan 470 menit di antaranya apnea — prevalensi 96,11%. Dataset ini tidak punya onset/durasi event maupun sleep stage.”

---

## 2. Membaca baris kedua: `ucddb / ucddb002`

| Kolom | Nilai | Artinya |
|---|---|---|
| `n_signals` | 14 | PSG 14 kanal |
| `signal_names` | `...ECG;SpO2;Sound;Flow;Sum;ribcage;abdo...` | Ada ECG, Flow, ribcage, abdo |
| `fs_ecg` | 128 | ECG 128 Hz |
| `fs_resp` | 8 | Flow 8 Hz |
| `duration_h` | 6.24 | 6,24 jam |
| `label_source` | `_stage.txt + _respevt.txt` | Stage + event list |
| `label_unit_s` | 30 | Satu label = epoch 30 detik |
| `n_label_units` | 748 | 748 epoch |
| `n_apnea_units` | 182 | 182 epoch beririsan dengan event |
| `apnea_prevalence` | 0.2433 | 182/748 = 24,33% |
| `n_events` | 124 | 124 event |
| `n_apnea_events` | 124 | Semuanya HYP/APNEA |
| `event_types` | `{"APNEA-M":2,"APNEA-O":7,"HYP-C":50,"HYP-M":3,"HYP-O":62}` | Rincian tipe |
| `ahi` | 23 | AHI klinis dari `SubjectDetails.xls` |
| `sleep_efficiency` | 84 | Sleep efficiency klinis (%) |
| `sleep_fraction` | 0.8382 | 83,82% epoch dinilai tidur |

**Kalimat satu baris:**
> “`ucddb002` adalah PSG 14 kanal selama 6,24 jam; ECG 128 Hz, Flow 8 Hz; 748 epoch 30 detik dengan 182 epoch apnea (24,33%); 124 event, AHI 23, sleep efficiency 84%.”

---

## 3. Hubungan CSV → JSON (buktikan sendiri)

JSON tidak ajaib; isinya penjumlahan kolom CSV.

**Prevalensi per dataset:**
- Apnea-ECG: `1609 / 3955 = 0.4068`
- slpdb: `3340 / 10197 = 0.3275`
- UCDDB: `4412 / 20789 = 0.2122`
- Overall: `9361 / 34941 = 0.2679`

**Total jam:**
- `66.1 + 85.62 + 173.36 = 325.08`

**Total record:**
- `8 + 18 + 25 = 51`

**Total subject:**
- `8 + 16 + 25 = 49`

**Prevalensi total:**
- `9361 / 34941 = 0.2679`

Kalau dosen tanya “angka 26,79% dari mana?”, jawab: **jumlah `n_apnea_units` dibagi jumlah `n_label_units` seluruh CSV.**

---

## 4. Membedakan `n_events` vs `n_apnea_events`

Contoh `slpdb / slp01a`:

- `n_events = 252` → semua token event: `H=23`, `HA=15`, `L=64`, `LA=150`.
- `n_apnea_events = 38` → hanya `H + HA = 23 + 15`.
- `L/LA` = leg movement, bukan apnea.
- `n_apnea_units = 37` → 37 epoch punya minimal satu event apnea. Bedanya dengan 38 karena satu epoch bisa punya lebih dari satu event.

Jadi:
- `n_events` = semua anotasi.
- `n_apnea_events` = hanya apnea/hypopnea.
- `n_apnea_units` = jumlah epoch/menit yang positif.

---

## 5. Kalau ditanya lebih dalam: jawaban cepat

| Pertanyaan | Jawaban |
|---|---|
| Kenapa 51 record tapi 49 subject? | `slp01a/b` dan `slp02a/b` subject yang sama |
| Kenapa Apnea-ECG cuma 8? | Hanya `*er` yang punya respirasi |
| Kenapa `n_events` Apnea-ECG kosong? | Label hanya A/N per menit, tanpa onset/durasi |
| Kenapa prevalensi UCDDB 21% tapi AHI 24? | Prevalensi = fraksi epoch; AHI = event/jam tidur |
| Apa itu `label_unit_s`? | Panjang satu unit label: 60 s atau 30 s |
| Apa itu `n_apnea_units`? | Jumlah unit label yang positif apnea |
| Apa itu `sleep_fraction`? | Fraksi epoch yang dinilai tidur |
| Apa beda `sleep_fraction` dan `sleep_efficiency`? | Fraksi epoch vs metrik klinis UCDDB |
| Kenapa `warnings` kosong? | Semua record lolos cek ECG + respirasi + label + durasi |
| Apakah prevalensi sama dengan AHI? | Tidak; beda satuan dan definisi |

---

## 6. Latihan 5 menit

1. Buka `outputs/dataset_inventory.csv` di Excel.
2. Filter `dataset = ucddb`.
3. Sort `apnea_prevalence` dari besar ke kecil.
   - Tertinggi: `ucddb025` = 0.8298 (83%).
   - Terendah: `ucddb018` = 0.0182 (1,8%).
4. Pilih satu baris, cek `n_apnea_units / n_label_units` dengan kalkulator; harus sama dengan `apnea_prevalence`.
5. Jumlahkan `n_label_units` untuk `ucddb`; hasilnya harus 20.789, sama dengan JSON.

Kalau lima langkah ini bisa kamu lakukan, kamu sudah paham cara bacanya.

---

## 7. Kalau benar-benar tidak tahu saat ditanya

Jangan mengarang. Jawab:

> “Saya tidak hafal angkanya, tapi kolom itu artinya [baca dari data dictionary], dan sumbernya dari [file mentah]. Saya bisa cek langsung di CSV/script.”

Yang penting kamu tahu **asal-usul** dan **cara verifikasi**, bukan menghafal 51 baris.
