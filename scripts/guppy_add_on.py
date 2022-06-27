# %%%

import seaborn as sns
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import h5py
import os
import numpy as np
import pandas as pd
import re


"""
general workflow for post processing is as follows::

retrieve data from guppy outputs:

    can be timestamps, FP data for plotting, static point data

for unprocessed timestamps, such as behavior we need to:
        do event aligning,

"""


"""
TODO class for subject data
TODO visualization module 
"""


class Subject:
    """
    The Subject class represents an experimental animal, it allows you to call any of the files from from the guppy output
    """

    def __init__(self, output_path):
        self.output_path = output_path

    def get_output_path(self):
        # returns the path to the output folder
        return self.output_path

    def get_raw_data(self):
        # reads storesnames file and reads all raw data files, returns as dictionary {storename: array}
        stores_list_file = os.path.join(self.output_path, 'storesList.csv')
        stores_list = pd.read_csv(stores_list_file).columns.tolist()
        raw_data = {}
        for i in os.listdir(self.output_path):
            for s in stores_list:
                if re.search(s, i):
                    read_h5 = h5py.File(os.path.join(path, i), 'r').get('data')
                    raw_data[s] = np.array(read_h5)
        return raw_data

    def get_PSTH_data(self, data_type='z_score'):
        regex = f"[^Uncorrected]_{data_type}_(.*?).h5"
        all_files = [i for i in os.listdir(
            self.output_path) if re.search(regex, i)]
        psth_files = []
        for i in all_files:
            if re.search("AUC", i):
                pass
            elif re.search('freq', i):
                pass
            else:
                psth_files.append(i)
        psth_data = {}
        for p in psth_files:
            key = p.split(f'_{data_type}_')[0]
            value = pd.read_hdf(os.path.join(self.output_path, p))
            psth_data[key] = value
        return psth_data


"""
Visualization module
"""


def plot_psth(data, x='timestamps', y='mean',  color=None, linewidth=None, alpha=0.2, ax=None):
    # function for line plots
    # if smooth = True:
    #     data = data.rolling(window=10).mean()
    # else:
    #     data = data
    sns.lineplot(data=data,
                 x=x,
                 y=y,
                 color=color,
                 linewidth=linewidth,
                 ax=ax
                 )
    if ax == None:  # will use plt object if axes is not declared
        plt.fill_between(data=data, x=data[x], y1=(data[y] + data['err']),
                         y2=(data[y] - data['err']), color=color, alpha=alpha, linewidth=0)
        sns.despine(top=True, right=True, ax=ax)
        ax.set_xlim(left=-5, right=15)  # x axis range
        ax.set_ylim(bottom=-0.5, top=2.5)  # y axis range
        ax.tick_params(axis="both", which="major",
                       labelsize=8)  # set tick label size
        ax.set_xlabel("Time(sec)", fontsize=8)  # x axis label
        ax.set_ylabel("z-score", fontsize=8)  # y axis label
        ax.set_box_aspect(1)
    else:  # will use axes object if axes is  declared
        ax.fill_between(data=data, x=data[x], y1=(data[y] + data['err']),
                        y2=(data[y] - data['err']), color=color, alpha=alpha, linewidth=0)
        sns.despine(top=True, right=True, ax=ax)
        ax.set_xlim(left=-5, right=15)  # x axis range
        ax.set_ylim(bottom=-0.5, top=2.5)  # y axis range
        ax.tick_params(axis="both", which="major",
                       labelsize=8)  # set tick label size
        ax.set_xlabel("Time(sec)", fontsize=8)  # x axis label
        ax.set_ylabel("z-score", fontsize=8)  # y axis label
        ax.set_box_aspect(1)


# %%
group = Subject(
    r'R:\Mike\LHA_dopamine\LH_NAC_Headfix_FP\Photometry\Pav Training\2_color_pav\Sucrose_to_sucralose\Training\Day_9\average')
group.get_output_path()
group_data = group.get_PSTH_data()
# %%
fig, axs = plt.subplots(nrows=2, ncols=1, figsize=(10, 10))
plot_psth(data=group_data['cue_GCAMPNAC'],
          color='black', linewidth=1, alpha=0.2, ax=axs[0])
plt.show()

# %%
sns.lineplot(data=group_data['cue_GCAMPNAC'], x='timestamps', y='mean')
plt.fill_between(data=group_data['cue_GCAMPNAC'], x=group_data['cue_GCAMPNAC']['timestamps'],
                 y1=(group_data['cue_GCAMPNAC']['mean'] +
                     group_data['cue_GCAMPNAC']['err']),
                 y2=(group_data['cue_GCAMPNAC']['mean'] -
                     group_data['cue_GCAMPNAC']['err']),
                 color='black', linewidth=0)

# %%


class Behavior:
    pass


mouse_1 = Subject(r'R:\Mike\LHA_dopamine\LH_NAC_Headfix_FP\Photometry\Pav Training\2_color_pav\Sucrose_to_sucralose\Training\Day_9\512581-220531-110552\512581-220531-110552_output_1')


# %%
for i in os.listdir(path):
    if re.search('AUC', i):
        print(i)
# %%


# %%
