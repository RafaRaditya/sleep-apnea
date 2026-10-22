# EDA Report — Sleep Apnea Multi-Dataset

## Dataset summary

```
{
  "apnea-ecg": {
    "n_records": 8,
    "n_subjects": 8,
    "total_duration_h": 66.1,
    "n_label_units": 3955,
    "n_apnea_units": 1609,
    "apnea_prevalence": 0.4068,
    "n_events": null,
    "n_apnea_events": null,
    "mean_ahi": null,
    "mean_sleep_fraction": null
  },
  "slpdb": {
    "n_records": 18,
    "n_subjects": 16,
    "total_duration_h": 85.62,
    "n_label_units": 10197,
    "n_apnea_units": 3340,
    "apnea_prevalence": 0.3275,
    "n_events": 4061,
    "n_apnea_events": 3496,
    "mean_ahi": null,
    "mean_sleep_fraction": 0.6875
  },
  "ucddb": {
    "n_records": 25,
    "n_subjects": 25,
    "total_duration_h": 173.36,
    "n_label_units": 20789,
    "n_apnea_units": 4412,
    "apnea_prevalence": 0.2122,
    "n_events": 3337,
    "n_apnea_events": 3318,
    "mean_ahi": 24.24,
    "mean_sleep_fraction": 0.7711
  },
  "overall": {
    "n_records": 51,
    "n_subjects": 49,
    "total_duration_h": 325.08,
    "n_label_units": 34941,
    "n_apnea_units": 9361,
    "apnea_prevalence": 0.2679
  }
}
```

## Window estimates

```
  dataset  n_records   n30  pos30   n10  pos10  prev30  prev10
apnea-ecg          8  7910   3218 23730   9654  0.4068  0.4068
    slpdb         18 10197   3340 30591  10020  0.3275  0.3275
    ucddb         25 20789   4412 62367   8982  0.2122  0.1440
```

## UCDDB: 30-s label propagated to 10-s windows

- 10-s windows: 62367
- Agreement with event-level 10-s labels: 0.9318
- Propagated positive but event negative: 4254
- Propagated negative but event positive: 0

## Sleep vs wake apnea prevalence

```
dataset  n_records  mean_sleep_fraction  mean_prev_sleep  mean_prev_wake
  slpdb         18                0.687            0.432           0.122
  ucddb         25                0.771            0.252           0.088
```

## UCDDB event duration

```
count    3318.0
mean       18.4
std         6.3
min         8.0
10%        12.0
50%        17.0
90%        27.0
max        65.0
```

Fraction of events shorter than 10 s: 0.006

## Implications

- **Window 10 s:** only UCDDB has event onsets and durations. For slpdb and Apnea-ECG, 10-s labels are inherited from 30-s/60-s units.
- **Wake/sleep:** slpdb and UCDDB have sleep stages; Apnea-ECG does not. Excluding wake changes the label distribution and cannot be applied consistently to Apnea-ECG.
- **Subject split:** slp01a/b and slp02a/b share subjects; group by subject, not record.