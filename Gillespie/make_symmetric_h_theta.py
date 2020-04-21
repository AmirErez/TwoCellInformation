#!/usr/bin/env python
#! Produce input script for simulating producer-consumer pairs
import pandas as pd
import numpy as np
from collections import defaultdict

outfile='symmetric_h_theta.txt'

ret = defaultdict(lambda: [])

# log10theta = np.linspace(-6, 0, 61);
# thetas = np.power(10, log10theta)
# thetas = np.linspace(-0.3, 0.3, 21)
hs = np.linspace(-0.3, 0.3, 61)
thetas = np.linspace(-0.3, 0.3, 61)

for hh in hs:
   for tt in thetas:
       ret['hx'].append(hh)
       ret['hy'].append(hh)
       ret['theta_x'].append(tt)
       ret['theta_y'].append(tt)

df = pd.DataFrame(ret)
df['log10nc_x'] = 3.4771
df['log10nc_y'] = 3.4771

df = df[['theta_x', 'theta_y', 'hx', 'hy', 'log10nc_x', 'log10nc_y']]
df.to_csv(outfile, sep=',', header=None, index=None, float_format='%g')
print('Finished successfully writing to {}'.format(outfile))
