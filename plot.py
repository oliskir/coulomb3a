import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


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


dist = 15


energies = [2000, 3000, 4000, 5000, 6000]
colors = ['C0', 'C1', 'C2', 'C3', 'C4']


arrows = []
labels = []


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
            lab = r'$E_{{2\alpha}}={0:.1f}$ MeV'.format(e/1000.)
            arr = plt.arrow(x1, y1, dx, dy, head_width=0.02, head_length=0.02, color=c, length_includes_head=True, label=lab)
            arrows.append(arr)
            labels.append(lab)
        else:
            plt.arrow(x1, y1, dx, dy, head_width=0.02, head_length=0.02, color=c, length_includes_head=True)

plt.xlim(-0.05,1.05)
plt.ylim(-0.05,1.05)

plt.plot([0, 0], [0, 1], color='black', linewidth=1, linestyle='dotted')

x = [0, np.cos(30*np.pi/180)]
y = [0, np.sin(30*np.pi/180)]
plt.plot(x, y, color='black', linewidth=1, linestyle='dotted')

q = (30 + 60 * np.linspace(0, 1, 100)) * np.pi/180
x = np.cos(q)
y = np.sin(q)
plt.plot(x, y, color='black', linewidth=1, linestyle='dotted')

plt.xlabel(r'$x = \sqrt{3} (E_2 - E_3) / E_{tot}$')
plt.ylabel(r'$y = (2 E_1 - E_2 - E_3) / E_{tot}$')

plt.title(r'$d_{{1}}={0:.0f}$ fm, $d_{{23}}=4.5$ fm, $E_{{3\alpha}}=9.35$ MeV'.format(dist))


for c,lab in zip(colors, labels):
    plt.scatter(np.nan, np.nan, c=c, marker=r'$\rightarrow$', s=100, label=lab)

plt.legend(prop={'size': 8})

plt.show()

