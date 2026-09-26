"""Download the datasets for the DK chapter 8 illustrations into data/.

The Rdatasets files are not committed: they come from R packages via the
Rdatasets mirror, whose redistribution terms are not definitive (see
data/README.md). The FRED yields are public domain and committed; the script
rebuilds them. Run from the repository root:

    .venv/bin/python data/fetch_dk_data.py
"""

import pathlib
import urllib.request

BASE = "https://raw.githubusercontent.com/vincentarelbundock/Rdatasets/master/csv"

# (local file, Rdatasets path, DK section)
DATASETS = [
    ("seatbelts.csv", "datasets/Seatbelts.csv", "8.2-8.3"),
    ("wwwusage.csv", "datasets/WWWusage.csv", "8.4"),
    ("mcycle.csv", "MASS/mcycle.csv", "8.5"),
]


# FRED H.15 monthly constant-maturity Treasury yields, as a public-domain
# substitute for DK 8.6's Diebold-Li (CRSP) data: (series, maturity in months)
FRED_URL = "https://fred.stlouisfed.org/graph/fredgraph.csv?id="
YIELDS = [("GS3M", 3), ("GS6M", 6), ("GS1", 12), ("GS2", 24), ("GS3", 36),
          ("GS5", 60), ("GS7", 84), ("GS10", 120)]
YIELDS_START, YIELDS_END = "1985-01-01", "2000-12-01"


def fetch_yields(dest):
    ids = ",".join(s for s, _ in YIELDS)
    with urllib.request.urlopen(FRED_URL + ids) as resp:
        lines = resp.read().decode().splitlines()
    rows = [l for l in lines[1:] if YIELDS_START <= l.split(",")[0] <= YIELDS_END]
    header = "date," + ",".join(f"m{m}" for _, m in YIELDS)
    dest.write_text("\n".join([header] + rows) + "\n")


def main():
    here = pathlib.Path(__file__).resolve().parent
    for name, path, section in DATASETS:
        dest = here / name
        if dest.exists():
            print(f"{name}: already present")
            continue
        with urllib.request.urlopen(f"{BASE}/{path}") as resp:
            dest.write_bytes(resp.read())
        print(f"{name}: downloaded (DK {section})")
    dest = here / "us_yields.csv"
    if dest.exists():
        print("us_yields.csv: already present")
    else:
        fetch_yields(dest)
        print("us_yields.csv: downloaded (DK 8.6 substitute)")


if __name__ == "__main__":
    main()
