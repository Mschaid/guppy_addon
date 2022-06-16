# %%
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.gridspec import GridSpec  # for subplots
import seaborn as sns
import pavlovian_functions as pv
import string

""" 
This script is used to create a dataframe of all the data from the pavlovian task.
all_data_path is the directory in which all data is stored.
new_folder is the name of the new folder to create. This is where all aggreated data will be saved.
"""
all_data_path = r'R:\Mike\LHA_dopamine\LH_NAC_Headfix_FP\Photometry\Pav Training\2_color_pav\Sucrose_to_sucralose\Training'
new_folder = 'aggregated_data'
pv.create_new_dir(all_data_path, new_folder)
path_to_save = os.path.join(all_data_path, new_folder)

# searches and collects all tidy data from day analysis
files = glob.glob(f"{all_data_path}/**/**/*tidy*.feather")

# splits tidy data into behavior. fiberphotometry, and AUC lists
all_behavior_files = [f for f in files if 'behavior' in f]
all_fp_files = [f for f in files if 'tidy_photometry' in f]
all_AUC_files = [f for f in files if 'AUC' in f]


"""
COMPILE, AGGREAGTE AND SAVE BEHAVIOR DATA
"""
behavior_df = pv.compile_data(all_behavior_files)
# all behavior from all days
behavior_group_by = ['time_sec', 'recording', 'day']

behavior_agg_dict = {'avg_frequency': ['mean', 'std', 'sem']}
# final dataframe of behavior
behavior_group_agg = (
    behavior_df
    .reset_index(drop=True)
    .groupby(by=behavior_group_by)
    .agg(behavior_agg_dict)
    .pipe(pv.flatten_df)
)
# save as feather file
behavior_group_agg.to_feather(
    f"{all_data_path}/{new_folder}/behavior_group_agg.feather")

"""
COMPILE, AGGREAGTE AND SAVE FIBERPHOTOMETRY DATA
"""

fp_df = pv.compile_data(all_fp_files)


def agg_and_clean_fp(df):
    """ function to agg and clean group fp data"""
    # left in script for piping purposes

    group_by = ['day', 'time_stamps', 'sensor', 'region', 'event']
    agg_dict = {'z-score': ['mean', 'std', 'sem']}

    df = (
        df
        .reset_index(drop=True)
        .groupby(by=group_by)
        .agg(agg_dict)
        .pipe(pv.flatten_df))
    return df


# select good nac mice, group and aggregate
nac_mice = ['512584', '514958']
final_nac_fp = (fp_df[fp_df.mouse.isin(nac_mice)]  # dropped missing nac probe mice
                .pipe(agg_and_clean_fp)
                )
final_nac_fp.to_feather(f"{all_data_path}/{new_folder}/final_nac_fp.feather")
final_nac_fp = final_nac_fp.iloc[::5, :]
# select good lha mice, group and aggregate
lha_mice = ['512581', '514959']
final_lha_fp = (fp_df[fp_df.mouse.isin(lha_mice)]  # dropped missing nac probe mice
                #  .pipe(drop_mice)
                .pipe(agg_and_clean_fp)
                )

final_lha_fp.to_feather(f"{all_data_path}/{new_folder}/final_lha_fp.feather")
final_lha_fp = final_lha_fp.iloc[::5, :]


"""
COMPILE, AGGREAGTE AND SAVE FIBERPHOTOMETRY AUC DATA
"""
AUC_df = pv.compile_data(all_AUC_files)
AUC_df.reset_index(drop=True).to_feather(
    f"{all_data_path}/{new_folder}/AUC_df.feather")

"""
PLOT DAYS
"""
days = ['3', '6', '9']
lha_gcamp = final_lha_fp.query(
    "region == 'LHA' & sensor == 'GCAMP' & event == 'cue' & day == @days")

nac_gcamp = final_nac_fp.query(
    "region=='NAC' & sensor=='GCAMP' & event=='cue' & day==@days")

lha_rda = final_lha_fp.query(
    "region=='LHA' & sensor=='RDA' & event=='cue' & day==@days")
nac_rda = final_nac_fp.query(
    "region=='NAC' & sensor=='RDA' & event=='cue' & day==@days")


"""
Plotting functions 
"""


def plot_day(df, colors, ax):
    sns.lineplot(data=df,
                 x='time_stamps',
                 y='z-score_mean',
                 palette=colors,
                 hue='day',
                 linewidth=0.5,
                 ax=ax).legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                               fontsize=8, borderaxespad=0, frameon=False, title='Day', title_fontsize=8)


def plot_day_err(df, day, color, ax):
    err = df.query("day==@day")
    ax.fill_between(data=err, x=err['time_stamps'], y1=(err['z-score_mean'] + err['z-score_sem']),
                    y2=(err['z-score_mean'] - err['z-score_sem']), color=color, alpha=0.5, linewidth=0)


def plot_behave_err(df, day, color, ax):
    err = df.query("day==@day")
    ax.fill_between(data=err, x=err['time_sec'], y1=(err['avg_frequency_mean'] + err['avg_frequency_sem']),
                    y2=(err['avg_frequency_mean'] - err['avg_frequency_sem']), color=color, alpha=0.5, linewidth=0)


def format_fp(ax):
    ax.set_xlim(left=-5, right=15)  # x axis range
    # axs.set_ylim(bottom=-0.5, top=2.5)  # y axis range
    ax.tick_params(axis="both", which="major",
                   labelsize=8)  # set tick label size
    ax.set_xlabel("Time(sec)", fontsize=8)  # x axis label
    ax.set_ylabel("z-score", fontsize=8)  # y axis label
    ax.set_box_aspect(1)
    sns.despine(top=True, right=True, ax=ax)


def format_behavior(ax, y_label):
    ax.set_xlim(left=-5, right=15)  # x axis range
    # axs.set_ylim(bottom=-0.5, top=2.5)  # y axis range
    ax.tick_params(axis="both", which="major",
                   labelsize=8)  # set tick label size
    ax.set_xlabel("Time(sec)", fontsize=8)  # x axis label
    ax.set_ylabel(y_label, fontsize=8)  # y axis label
    ax.set_box_aspect(1)
    sns.despine(top=True, right=True, ax=ax)


# %%
fig = plt.figure(figsize=(8, 8))
days = ['3', '6', '9']
gcamp_colors = ['#c8ffc3', '#73ff68', '#17e006']
rda_colors = ['#ffafaf', '#ff6f6f', '#f11717']
behavior_colors = ['#BEBEBE', '#707070', '#000000']

gs = GridSpec(3, 4, figure=fig)

"""
GCAMP
"""
ax1 = fig.add_subplot(gs[0, 0])
plot_day(lha_gcamp, gcamp_colors, ax1)
ax1.get_legend().remove()
for day, color in zip(days, gcamp_colors):
    plot_day_err(lha_gcamp, day, color, ax1)
format_fp(ax1)
ax1.set_title('LHA DA terminal activty', fontsize=8)

ax2 = fig.add_subplot(gs[0, 1])
plot_day(nac_gcamp, gcamp_colors, ax2)
for day, color in zip(days, gcamp_colors):
    plot_day_err(nac_gcamp, day, color, ax2)
format_fp(ax2)
ax2.set_title('NAc DA terminal activty', fontsize=8)

"""
rDlight
"""
ax3 = fig.add_subplot(gs[1, 0])
plot_day(lha_rda, rda_colors, ax3)
ax3.get_legend().remove()
for day, color in zip(days, rda_colors):
    plot_day_err(lha_rda, day, color, ax3)
format_fp(ax3)
ax3.set_title('Dopamine in LHA', fontsize=8)

ax4 = fig.add_subplot(gs[1, 1])
plot_day(nac_rda, rda_colors, ax4)
for day, color in zip(days, rda_colors):
    plot_day_err(nac_rda, day, color, ax4)
format_fp(ax4)
ax4.set_title('Dopamine in NAc', fontsize=8)

"""
behavior
"""
ax5 = fig.add_subplot(gs[2, 0])
licks = behavior_group_agg.query("recording=='lick' & day==@days")

sns.lineplot(data=licks, x='time_sec', y='avg_frequency_mean', linewidth=0.5,
             hue='day', palette=behavior_colors)
ax5.get_legend().remove()
for day, color in zip(days, behavior_colors):
    plot_behave_err(licks, day, color, ax5)
format_behavior(ax5, 'Lick Rate')
ax5.set_title('Lick Frequency', fontsize=8)


ax6 = fig.add_subplot(gs[2, 1])
encoder = behavior_group_agg.query("recording=='encoder' & day==@days")

sns.lineplot(data=encoder, x='time_sec', y='avg_frequency_mean', linewidth=0.5,
             hue='day', palette=behavior_colors).legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                                                        fontsize=8, borderaxespad=0, frameon=False, title='Day', title_fontsize=8)
for day, color in zip(days, behavior_colors):
    plot_behave_err(encoder, day, color, ax6)
format_behavior(ax6, 'Velocity')
ax6.set_title('Movement', fontsize=8)

axes = [ax1, ax2, ax3, ax4, ax5, ax6]
for ax in axes:
    ax.yaxis.set_label_coords(-0.35, 0.5)


gs.update(left=0.1, right=0.9, top=0.965, bottom=0.03, wspace=0.8, hspace=-0.4)

# draw box for cue

for ax in axes:
    pv.draw_cue_box(ax, alpha=0.25, color='lightgrey')

# label plots with letters
for n, ax in enumerate(axes):
    ax.text(-0.5, 1.2, string.ascii_uppercase[n], transform=ax.transAxes,
            size=14, weight='bold')

plt.rcParams['svg.fonttype'] = 'none'  # save text as text in svg
plt.savefig(path_to_save + '\\day_3_6_9.tiff',
            dpi=300,
            transparent=True)
plt.savefig(path_to_save + '\\day_3_6_9.svg',
            dpi=300,
            transparent=True)
plt.show()


# %%
