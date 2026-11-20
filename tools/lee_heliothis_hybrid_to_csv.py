"""Convert K. P. Lee's workbook for the Heliothis hybrid (H. subflexa female x
H. virescens male; Hybrid.xls, supplied by K. P. Lee) to the long format of
data/raw/lee_heliothis_larvae_raw.csv.

Sheets used:
  nochoice  block 'protein-carbohydrate' (columns 120-135): cumulative protein
            and carbohydrate eaten (mg) per larva from day 0 to days 1-7, on the
            five foods p35:c7 ... p7:c35 (12 larvae per food);
  choice    block at columns 14-17: cumulative protein and carbohydrate (mg)
            per larva, days 0-1 ... 0-4, choice between p35:c7 and p14:c28.
Protein and carbohydrate were computed by Lee from dry mass eaten and the food
composition (e.g. protein = 0.35 x dry mass eaten on p35:c7).

Usage:
    python3 tools/lee_heliothis_hybrid_to_csv.py Hybrid.xls data/raw/lee_heliothis_hybrid_raw.csv
(needs pandas and xlrd)
"""
import sys

import pandas as pd

SPECIES = "Heliothis subflexa x virescens"
FOODS = {1: "p35:c7", 2: "p28:c14", 3: "p21:c21", 4: "p14:c28", 5: "p7:c35"}
x = pd.read_excel(sys.argv[1], sheet_name=None, header=None)
rows = []

t = x["nochoice"]
cur = None
for i in range(2, t.shape[0]):
    tr = t.iat[i, 120]
    if isinstance(tr, (int, float)) and pd.notna(tr):
        cur = int(tr)
    rep = t.iat[i, 121]
    if cur is None or pd.isna(rep) or isinstance(rep, str):
        continue
    for k in range(7):
        rows.append(dict(treat=FOODS[cur], replicate=int(rep), species=SPECIES, experiment="Nochoice",
                         day=f"day0-{k + 1}", day2=k + 1,
                         protein=pd.to_numeric(t.iat[i, 122 + 2 * k], errors="coerce"),
                         carb=pd.to_numeric(t.iat[i, 123 + 2 * k], errors="coerce")))

c = x["choice"]
day = None
for i in range(1, c.shape[0]):
    lab = c.iat[i, 14]
    if isinstance(lab, str) and lab.startswith("day0-"):
        day = int(lab.split("-")[1])
    rep = c.iat[i, 15]
    if day is None or pd.isna(rep) or isinstance(rep, str):
        continue
    rows.append(dict(treat="p35:c7 vs p14:c28", replicate=int(rep), species=SPECIES, experiment="Choice",
                     day=f"day0-{day}", day2=day,
                     protein=pd.to_numeric(c.iat[i, 16], errors="coerce"),
                     carb=pd.to_numeric(c.iat[i, 17], errors="coerce")))

d = pd.DataFrame(rows)
d.to_csv(sys.argv[2], index=False)
print(f"wrote {len(d)} rows ({d.groupby('experiment').replicate.nunique().to_dict()} larvae per treatment) to {sys.argv[2]}")
