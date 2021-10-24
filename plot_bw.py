import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats as stats


bins=50


df1 = pd.read_csv('output/bw_r10fm_N5000.dat', header=None, names=['exc0','exc'])
df2 = pd.read_csv('output/bw_r10fm_N5000_ad.dat', header=None, names=['exc0','exc'])

x = df1['exc0']
f = stats.gaussian_kde(x)
n, x, _ = plt.hist(x, bins=bins, range=(1.5,6.0), histtype='step', density=True, color='white')
plt.plot(x, f(x), label='No Coulomb repulsion')

x = df1['exc']
f = stats.gaussian_kde(x)
n, x, _ = plt.hist(x, bins=bins, range=(1.5,6.0), histtype='step', density=True, color='white')
plt.plot(x, f(x), label=r'Isotropic; $d_1=10$ fm')

x = df2['exc']
f = stats.gaussian_kde(x)
n, x, _ = plt.hist(x, bins=bins, range=(1.5,6.0), histtype='step', density=True, color='white')
plt.plot(x, f(x), label=r'$1-\sin^2\theta$; $d_1=10$ fm')


plt.xlabel(r'$E_{x}$(Be) (MeV)')

plt.title(r'Effect of Coulomb repulsion on reconstruction of $^8$Be excitation energy')

plt.legend(prop={'size': 8})

plt.savefig('bw.png')

plt.show()

