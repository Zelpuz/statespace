# Data

| File | Committed | Series | Source | Used for |
|---|---|---|---|---|
| `nile.csv` | yes | Annual flow of the Nile at Aswan, 1871–1970 | Cobb (1978). Public domain, as in `statsmodels.datasets.nile` | DK ch. 2, tests |
| `seatbelts.csv` | no | Monthly UK car drivers and front and rear seat passengers killed or seriously injured, with km driven, real petrol price and the seat belt law, 1969–1984 | Harvey and Durbin (1986), R `datasets::Seatbelts` | DK 8.2–8.3 |
| `wwwusage.csv` | no | Users logged on to an Internet server each minute, over 100 minutes | Makridakis, Wheelwright and Hyndman (1998), R `datasets::WWWusage` | DK 8.4 |
| `us_yields.csv` | yes | Monthly US Treasury constant-maturity yields, 1985–2000, at 3, 6, 12, 24, 36, 60, 84 and 120 months (FRED `GS3M`, `GS6M`, `GS1`, `GS2`, `GS3`, `GS5`, `GS7`, `GS10`) | Board of Governors of the Federal Reserve System, H.15 Selected Interest Rates, retrieved from FRED, Federal Reserve Bank of St. Louis. Public domain, citation requested | DK 8.6 (substitute) |
| `mcycle.csv` | no | Head acceleration against time after a simulated motorcycle impact | Silverman (1985), R `MASS::mcycle` | DK 8.5 |

The uncommitted files come from the [Rdatasets](https://github.com/vincentarelbundock/Rdatasets)
mirror. Its maintainer believes these datasets are free to redistribute but couldn't
confirm their terms, so we download them instead of committing them:

    .venv/bin/python data/fetch_dk_data.py

DK 8.6 uses the Diebold and Li (2006) yields: 17 maturities from the CRSP government
bond files, which are licensed data. `us_yields.csv` covers the same months with the
8 public-domain H.15 maturities instead. The dynamic factor example therefore shows
the method but doesn't reproduce DK's estimates.

Full citations for every dataset, including the one FRED requests, are in
[`docs/references.md`](../docs/references.md#data).
