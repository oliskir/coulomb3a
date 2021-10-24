import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams['text.usetex'] = True
plt.rcParams['font.size'] = 11
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'


def xy_dalitz(row):
    ei = row['a1']
    ej = row['a2']
    ek = row['a3']
    e = np.sort(np.array([ei, ej, ek]))
    e[:] = e[::-1]
    etot = np.sum(e)
    x = np.sqrt(3) * (e[1] - e[2]) / etot
    y = (2 * e[0] - e[1] - e[2]) / etot
    return x, y


dist = 10
show_leg = True
figname = f'coul3a_{dist}fm_pretty'


energies = [2000, 3000, 4000, 5000, 6000]
colors = ['C0', 'C1', 'C2', 'C3', 'C4']


arrows = []
labels = []


fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(6, 5))


for e,c in zip(energies, colors):

    df1 = pd.read_csv(f'output/energies_e{e}keV_r1000fm.dat', header=None, names=['angle','a1','a2','a3'])
    df2 = pd.read_csv(f'output/energies_e{e}keV_r{dist}fm.dat', header=None, names=['angle','a1','a2','a3'])

    df1 = df1[(df1.angle>=0) & (df1.angle<=90)]
    df2 = df2[(df2.angle>=0) & (df2.angle<=90)]

    for i in range(len(df1)):

        x1, y1 = xy_dalitz(df1.iloc[i])
        x2, y2 = xy_dalitz(df2.iloc[i])
        dx = x2 - x1
        dy = y2 - y1
    
        if i==0:
            lab = r'$E_{{23}}={0:.1f}$ MeV'.format(e/1000.)
            arr = ax.arrow(x1, y1, dx, dy, head_width=0.02, head_length=0.02, color=c, length_includes_head=True, label=lab)
            arrows.append(arr)
            labels.append(lab)
        else:
            ax.arrow(x1, y1, dx, dy, head_width=0.02, head_length=0.02, color=c, length_includes_head=True)

ax.set_xlim(-0.05,1.05)
ax.set_ylim(-0.05,1.05)

ax.plot([0, 0], [0, 1], color='black', linewidth=1, linestyle='dotted')

x = [0, np.cos(30*np.pi/180)]
y = [0, np.sin(30*np.pi/180)]
ax.plot(x, y, color='black', linewidth=1, linestyle='dotted')

q = (30 + 60 * np.linspace(0, 1, 100)) * np.pi/180
x = np.cos(q)
y = np.sin(q)
ax.plot(x, y, color='black', linewidth=1, linestyle='dotted')

ax.set_xlabel(r'$x = \sqrt{3} (E_2 - E_3) / \sum_{i} E_{i}$')
ax.set_ylabel(r'$y = (2 E_1 - E_2 - E_3) / \sum_{i} E_{i}$')

ax.text(0.7, 0.85, r'$d_{{1}}={0:.0f}$ fm'.format(dist) + '\n' + r'$d_{{23}}=4.5$ fm' + '\n' + r'$\sum_{i} E_{i}=9.35$ MeV')


for c,lab in zip(colors, labels):
    ax.scatter(np.nan, np.nan, c=c, marker=r'$\rightarrow$', s=100, label=lab)

if show_leg: ax.legend(loc='lower right', prop={'size': 10}, frameon=False)

ax.tick_params(top=True, right=True, which='both')

plt.savefig(figname+'.pdf')
plt.savefig(figname+'.png', dpi=600)

plt.show()

plt.close(fig)

