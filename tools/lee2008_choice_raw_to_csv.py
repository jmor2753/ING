"""Convert K. P. Lee's raw choice-experiment workbook (Lee et al. 2008 PNAS,
Drosophila; supplied by K. P. Lee) to a long CSV.

The workbook gives, for each fly, CUMULATIVE protein and carbohydrate intake
(ug) from day 0 to days 3, 6, ..., 24, for nine yeast x sucrose treatments
(5 blocks x 5 replicates). Two cells are marked "do not use" by the author;
they are kept with do_not_use = 1 and excluded downstream.

Usage:
    python3 tools/lee2008_choice_raw_to_csv.py Choice_data_Lee_et_al_2008_PNAS.xlsx data/lee2008_choice_raw.csv
"""
import sys

import pandas as pd

TREAT = {1: (180, 180), 2: (180, 90), 3: (180, 45), 4: (90, 180), 5: (90, 90),
         6: (90, 45), 7: (45, 180), 8: (45, 90), 9: (45, 45)}   # (yeast, sucrose) g/L
DAYS = [3, 6, 9, 12, 15, 18, 21, 24]

t = pd.read_excel(sys.argv[1], header=None)
rows = []
for i in range(4, t.shape[0]):
    b, tr, rep = t.iat[i, 4], t.iat[i, 5], t.iat[i, 6]
    if pd.isna(b):
        continue
    y, s = TREAT[int(tr)]
    for k, day in enumerate(DAYS):
        P, C = t.iat[i, 7 + 2 * k], t.iat[i, 8 + 2 * k]
        flag = int(isinstance(P, str) or isinstance(C, str))
        rows.append(dict(block=int(b), treatment=int(tr), replicate=int(rep),
                         yeast_gL=y, sucrose_gL=s, day=day,
                         P_ug_cum=pd.to_numeric(P, errors="coerce"),
                         C_ug_cum=pd.to_numeric(C, errors="coerce"),
                         do_not_use=flag))
pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
print(f"wrote {len(rows)} rows to {sys.argv[2]}")
