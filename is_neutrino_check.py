import numpy as np
import uproot
#import pandas as pd

f = uproot.open('/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.flat.root')
events = uproot.open('/pnfs/icarus/scratch/users/mueller/systematics/sample_cv.flat.root:recTree')
#events = events.arrays(["rec.dlp_true.is_neutrino"])

for batch in events.iterate(step_size="100 kB", library="np"):
    for i in batch['rec.dlp_true.is_neutrino']:
        print(i)
