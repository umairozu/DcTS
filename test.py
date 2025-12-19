from new_DcTS_cli import *
from pathlib import Path



model = CassetteTapeDecay.from_rawdata_xlsx("RawData.xlsx")

for temp in [60, 65, 70]:
    frac_E_3w = model.remaining_fraction(temp, encapsulated=True, weeks=3)
    frac_D_3w = model.remaining_fraction(temp, encapsulated=False, weeks=3)
    print(temp, "E 3w:", frac_E_3w, "D 3w:", frac_D_3w)


#<--------------------------------------------------------->


