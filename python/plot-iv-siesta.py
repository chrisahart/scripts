import numpy as np
import matplotlib.pyplot as plt
from scripts.main import parameters as param

""" Plotting of SMEAGOL output _TRC.agr by filename"""

plotting_colors = ['r', 'g', 'b', 'm', 'grey', 'orange', 'cyan']

xlim = [-1, 1]  # Li and Au chain
labels = ['0', '-1', '-2', '-1.25', '-1.5', '-1.75']
folder = ['/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-0',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-1',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-2',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-1.5',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-1.25',
          '/Volumes/Storage/Data/Work/Postdoc/Work/calculations/transport/iv/siesta/iv-bottom-1.75']

data_iv = []
data_chr = []
data_chr_end = []
for i in range(len(folder)):
    data_iv.append(np.genfromtxt('{}/Liwire.CUR'.format(folder[i]), skip_header=0, skip_footer=0))
    data_chr.append(np.genfromtxt('{}/Liwire.CHR'.format(folder[i]), skip_header=0, skip_footer=0))
    data_chr_end.append(data_chr[i][-1, 2])

# IV curve
fig_plot_1, ax_plot_1 = plt.subplots()
for i in range(len(folder)):
    ax_plot_1.plot(data_iv[i][:, 0], data_iv[i][:, 1], '.-', color=plotting_colors[i], label=labels[i])
ax_plot_1.legend(frameon=False)
ax_plot_1.set_xlabel('Bias voltage / eV')
ax_plot_1.set_ylabel('Current / A')
fig_plot_1.tight_layout()

# Charge
fig_plot_2, ax_plot_2 = plt.subplots()
ax_plot_2.plot(np.array(labels, dtype=np.float32), data_chr_end, 'k.')
ax_plot_2.axhline(y=12, color='r', linestyle='-')
ax_plot_2.set_xlabel('HartreeLeadsBottom')
ax_plot_2.set_ylabel('Charge difference')
fig_plot_2.tight_layout()

if __name__ == "__main__":
    print('Finished.')
    plt.show()
