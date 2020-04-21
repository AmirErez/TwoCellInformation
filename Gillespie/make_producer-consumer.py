#!/usr/bin/env python
#! Produce input script for simulating producer-consumer pairs
import pandas as pd
import numpy as np
from collections import defaultdict

outfile='producer-consumer.txt'

ret = defaultdict(lambda: [])

log10theta = np.linspace(-6, 0, 61);
thetas = np.power(10, log10theta)

# for theta in np.linspace(0,0.2,41):
for theta in np.power(10, log10theta):
   for hx in np.linspace(0.0,0.3,31):
#   for hx in [0.0]:
       hy = -hx
       ret['hx'].append(hx)
       ret['hy'].append(-hx)
       ret['theta_x'].append(theta)
       ret['theta_y'].append(theta)

df = pd.DataFrame(ret)
df['log10nc_x'] = 3.5
df['log10nc_y'] = 3.5

df = df[['theta_x', 'theta_y', 'hx', 'hy', 'log10nc_x', 'log10nc_y']]
df.to_csv(outfile, sep=',', header=None, index=None, float_format='%g')
print('Finished successfully writing to {}'.format(outfile))
