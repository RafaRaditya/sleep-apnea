# Bukti Validasi Rekaman & Harmonisasi Awal

Dokumen ini memetakan setiap klaim ke file/angka/perintah yang bisa ditunjukkan ke dosen.

## Ringkasan klaim

| Klaim | Angka | Bukti utama |
|---|---|---|
| Total rekaman | 51 record, 49 subject, 325,08 jam | `outputs/dataset_inventory.csv`, `outputs/dataset_summary.json` |
| Apnea-ECG yang punya respirasi | 8 record (`a01er`–`c03er`) | `apnea-ecg/DESTINATION/*er.hea`, kolom `resp_channel` |
| slpdb | 18 record, 16 subject | `slpdb/DESTINATION/RECORDS`, inventory |
| UCDDB | 25 record | `ucddb/DESTINATION/RECORDS`, inventory |
| Resolusi label | Apnea-ECG 60 s; slpdb 30 s; UCDDB event onset+duration | kolom `label_unit_s`, `label_source`; `slpdb/README.st`; `ucddb/*_respevt.txt` |
| Window 10 detik | UCDDB 10-s prevalence 14,40% vs 30-s 21,22% | `outputs/eda_report.md`, `outputs/window_estimates_by_record.csv` |
| Propagasi 30→10 detik | agreement 93,18%; 4.254 false positive; 0 false negative | `outputs/eda_report.md` |
| Sleep vs wake | slpdb 43,2% vs 12,2%; UCDDB 25,2% vs 8,8% | `outputs/sleep_wake_stats.csv`, `outputs/eda_report.md` |
| Durasi event UCDDB | median 17 s; 90% ≤27 s; 0,6% <10 s | `outputs/ucddb_event_durations.csv`, `outputs/figures/ucddb_event_durations.png` |
| Contoh sinyal + label | ECG + respirasi, apnea vs normal | `outputs/figures/example_windows.png` |

## 1. Membuktikan 51 record / 49 subject

```powershell
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_inventory.csv'); print('rows=', len(df)); print(df.groupby('dataset').size()); print('subjects=', df['subject_id'].nunique()); print('hours=', round(df['duration_s'].sum()/3600, 2))"
```

Hasil yang harus muncul: `rows= 51`, `apnea-ecg 8`, `slpdb 18`, `ucddb 25`, `subjects= 49`, `hours= 325.08`.

## 2. Membuktikan Apnea-ECG hanya 8 record yang punya respirasi

```powershell
Get-Content apnea-ecg\DESTINATION\a01er.hea
Get-ChildItem apnea-ecg\DESTINATION\a01er.* | Select-Object Name, Length
```

Header `a01er.hea` menunjukkan 5 sinyal: `ECG`, `Resp C`, `Resp A`, `Resp N`, `SpO2`. File `a01er.dat` tidak ada; ECG ada di `a01.dat` dan respirasi di `a01r.dat`. Ini sebabnya loader membaca keduanya terpisah.

Cek di inventory:

```powershell
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_inventory.csv'); print(df[df.dataset=='apnea-ecg'][['record_id','signal_names','resp_channel']].to_string(index=False))"
```

## 3. Membuktikan slpdb punya respirasi dan label 30 detik

```powershell
Get-Content slpdb\DESTINATION\slp01a.hea
Get-Content slpdb\DESTINATION\README.st -TotalCount 25
```

`slp01a.hea` menunjukkan `ECG`, `BP`, `EEG (C4-A1)`, `Resp (sum)`. `README.st` mendefinisikan stage `W/1/2/3/4/R` dan event `H/HA/OA/X/CA/CAA/L/LA/A/MT`.

Cek resolusi label:

```powershell
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_inventory.csv'); print(df[df.dataset=='slpdb'][['record_id','label_unit_s','n_label_units','n_apnea_units']].to_string(index=False))"
```

## 4. Membuktikan UCDDB harus memakai `.rec`, bukan `_lifecard.edf`

```powershell
python -c "import pyedflib; f=pyedflib.EdfReader(r'ucddb\DESTINATION\ucddb002.rec'); print('REC:', f.getSignalLabels()); f.close()"
python -c "import mne; raw=mne.io.read_raw_edf(r'ucddb\DESTINATION\ucddb002_lifecard.edf', preload=False, verbose='ERROR'); print('LIFECARD channels:', raw.ch_names, 'n=', len(raw.ch_names))"
```

Hasil: `.rec` punya 14 kanal termasuk `ECG`, `Flow`, `ribcage`, `abdo`, `SpO2`; `_lifecard.edf` hanya 3 kanal ECG (`chan 1`–`chan 3`). Catatan: pyedflib menolak header `_lifecard.edf` karena format starttime-nya tidak standar; MNE bisa membacanya.

Cek label UCDDB:

```powershell
Get-Content ucddb\DESTINATION\ucddb002_respevt.txt -TotalCount 8
python -c "import pandas as pd; df=pd.read_csv('outputs/dataset_inventory.csv'); print(df[df.dataset=='ucddb'][['record_id','label_unit_s','n_events','n_apnea_events','ahi']].head(5).to_string(index=False))"
```

`_respevt.txt` menunjukkan kolom `Time`, `Type`, `Duration` — inilah yang membuat UCDDB satu-satunya dataset dengan label event-level 10 detik.

## 5. Membuktikan klaim window 10 detik

```powershell
Get-Content outputs\eda_report.md
```

Lihat bagian:

- **Window estimates**: UCDDB `prev30=0.2122`, `prev10=0.1440`.
- **UCDDB: 30-s label propagated to 10-s windows**: agreement `0.9318`, false positive `4254`, false negative `0`.
- **Sleep vs wake apnea prevalence**: slpdb `0.432` vs `0.122`; UCDDB `0.252` vs `0.088`.
- **UCDDB event duration**: median `17 s`; fraction <10 s `0.006`.

File pendukung: `outputs/window_estimates_by_record.csv` (satu baris per record: `n30`, `pos30`, `n10`, `pos10`, `prev30`, `prev10`).

## 6. Membuktikan loader benar-benar membaca data

```powershell
python code\dataset_inventory.py --self-test
```

Output yang harus muncul: tiga blok (`apnea-ecg/a01er`, `slpdb/slp01a`, `ucddb/ucddb002`) dengan shape sinyal, fs, jumlah label, dan contoh event. Ini membuktikan bukan data karangan: loader membaca file mentah WFDB/EDF.

## 7. Membuktikan kualitas sinyal

Buka `outputs/figures/example_windows.png`. Terlihat:
- Apnea-ECG punya artefak besar di awal rekaman (ECG ±10 mV).
- UCDDB Flow menyentuh batas ±20 L/min.
- slpdb relatif bersih.

Ini alasan preprocessing/clipping masuk sebagai kebutuhan, bukan sekadar ablation.

## 8. Membuktikan CSV/JSON tidak salah hitung

```powershell
python code\verify_inventory.py --checksums all
```

Script ini melakukan dua hal:
1. **Menghitung ulang** `n_signals`, `duration_s`, `n_label_units`, `n_apnea_units`, `n_events`, `n_apnea_events`, dan `apnea_prevalence` langsung dari file mentah WFDB/EDF, lalu membandingkannya dengan `dataset_inventory.csv`.
2. **Mengecek SHA256** file mentah yang dipakai terhadap `SHA256SUMS.txt` resmi PhysioNet.

Hasil yang sudah didapat:

```text
Records: 51 | PASS: 51 | FAIL: 0
JSON aggregate problems: 0
Checksum problems (all): 0
```

Output: `outputs/verification_report.md` dan `outputs/verification_results.csv` (status per record).

## Urutan demo 2 menit

1. Buka `outputs/dataset_inventory.csv` di Excel → filter `dataset`. Tunjukkan 51 baris dan kolom `resp_channel`, `label_unit_s`, `n_apnea_units`.
2. Buka `outputs/dataset_summary.json` → tunjukkan total 325,08 jam dan 49 subject.
3. Buka `outputs/eda_report.md` → tunjukkan tabel window dan bagian propagasi 30→10 detik.
4. Buka `outputs/figures/example_windows.png` → tunjukkan contoh apnea vs normal.
5. Jalankan `python code\verify_inventory.py --checksums all` → tunjukkan **51/51 PASS, 0 mismatch, 0 checksum problem**.
6. Jalankan `python code\dataset_inventory.py --self-test` → bukti loader membaca file mentah.
7. Kalau ditanya “sumbernya dari mana?”, buka `apnea-ecg\DESTINATION\a01er.hea`, `slpdb\DESTINATION\slp01a.hea`, dan `ucddb\DESTINATION\ucddb002.rec` (atau jalankan perintah pyedflib di atas).

## Kalau ditanya “ini baru baca dataset, bukan riset?”

Jawab: ini **tahap validasi rekaman**, dan pertanyaan riset pertama di proposal memang harmonisasi multi-dataset. Audit ini menghasilkan temuan yang mengubah desain eksperimen:

- Hanya 8 record Apnea-ECG yang layak multimodal.
- UCDDB harus memakai `.rec`; kalau tidak, dataset gagal validasi.
- Window 10 detik hanya bersih di UCDDB; di dataset lain labelnya warisan.
- Wake/sleep mengubah prevalensi 2–3× lipat.

Tanpa temuan ini, model bisa dilatih di atas label yang tidak sebanding antar-dataset. Tahap berikutnya: harmonisasi sinyal (resampling, filtering, normalisasi) + split per subject + ekstraksi window.
