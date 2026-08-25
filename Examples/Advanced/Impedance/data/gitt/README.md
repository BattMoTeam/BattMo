# BAK N21700CG GITT example data

These files support the BattMo GITT-to-ECM calibration example.

- `gitt_measurements.csv` contains measured GITT time-series data for a nominal 5 Ah BAK N21700CG cell.
- `gitt_ocv_soc.csv` contains the vetted directional OCV reference extracted by the legacy PyBOP workflow.
- `gitt_rc_parameters.csv` contains the vetted legacy per-pulse 1-RC fits.
- `gitt_ecm_parameters.csv` contains the vetted legacy interpolated 1-RC parameter tables.

The source data uses positive current for charge. BattMo's equivalent-circuit model uses positive current for discharge, so the example loader changes the sign exactly once.

The files were copied from the Task-5.4 repository at commit `3a4a04157fa312978a89a9c425c7ee3200398b96`. Redistribution in BattMo was authorized for this example. The calibration code does not use the reference parameters as fitting inputs; they are retained for regression comparison.

The measurement CSV is a stand-alone, lossless column selection from
`DigiBatt-BAK-5000-N21700CG-006-GITT-data.parquet`. It contains all 560779
measurement rows and the five source columns used by this example:
`Time_h`, `U_V`, `I_A`, `DataSet`, and `T1_C`. The source Parquet is not
required or read by BattMo. The conversion was performed separately with
DuckDB's `read_parquet` and `COPY ... (HEADER, DELIMITER ',')` commands.

## SHA-256

| File | SHA-256 |
| --- | --- |
| `gitt_measurements.csv` | `77774e983694d3dab009a809741a2c9bbd94b4f14c76134abac7cb38841b2bf8` |
| `gitt_ecm_parameters.csv` | `8a3d484bbb4b40c13959640efccd8c6a1fb602c87d657fd8a1002d8ceb74ee53` |
| `gitt_ocv_soc.csv` | `348be5cec4841b6e073d655ff410f19f87ad139eefe695fe19cec9de283c7e07` |
| `gitt_rc_parameters.csv` | `24c6ef52a43cb33a237a779d337bd20c2cfcc5249e4a799c13eff3d6bb2b175a` |

The SHA-256 of the source Parquet before conversion was
`39704b7a26977784dd105102b142f063cf064bcdf446eeb599ad69e117eca62d`.

## Cell metadata

- Cell: BAK N21700CG
- Nominal capacity: 5 Ah
- Nominal voltage: 3.6 V
- Charge limit: 4.2 V
- Discharge limit: 2.5 V
- Temperature column: `T1_C`
