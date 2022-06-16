
import numpy as np
import pandas as pd
import h5py
import os
import glob
import re
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import seaborn as sns


def create_new_dir(file_path, new_dir_ext):
    """ creates new directory from file path and extension"""
    new_directory = os.path.join(file_path, new_dir_ext)
    if not os.path.exists(new_directory):
        os.mkdir(new_directory)
        print('directory created')
    else:
        print('directroy already exists')


def find_file_path(file_path, file_basename):
    """ find full path when given filepath and basename of file"""
    for root, dirs, files in os.walk(file_path):
        if file_basename in files:
            return os.path.join(root, file_basename)


def read_timestamps(file_path):
    """ reads raw timesstamp files from hdf5 guppy foles"""
    ts_arr = h5py.File(file_path, "r").get('timestamps')
    return np.array(ts_arr)


def list_subdirs(directory):
    """ given directory, returns list of all sub dirs as full path"""
    return [os.path.join(directory, file) for file in os.listdir(directory)]


def align_events(df, event, event_align):
    """ function to create PSTH by aligning events to specific event, given a dataframe"""
    event_arr = df[event][np.where(df[event] != 0)[0]].dropna()
    event_align_arr = df[event_align].dropna().to_numpy()
    align_events_dict = {}
    for i in range(event_arr.shape[0]):
        arr = np.array(event_align_arr - event_arr[i])
        align_events_dict[i] = arr
        new_df = pd.DataFrame(dict([(k, pd.Series(v))
                              for k, v in align_events_dict.items()]))
    return new_df


def calc_frequency(filepath):
    """ calculate frequency of events around specifc EPOCH"""
    # time is in seconds
    df = pd.read_feather(filepath)  # read csv
    arr = df[(df > -10) & (df < 20)].to_numpy()  #
    arr = arr[np.logical_not(np.isnan(arr))]
    freq = np.histogram(arr, bins=155)[0] / len(df.columns)
    freq_convert = freq * 5
    return freq_convert


def flatten_df(df):
    """ flattens column to tidy data after groupby and aggreation"""
    df.columns = df.columns = ['_'.join(col) for col in df.columns.values]
    df = df.reset_index()
    return df


def create_behavior_dataframe(recording: str, filepath_list):  # recording: str,
    """ function to calculate mean frequency of event and convert to dataframe"""

    freq_dict = {}
    for f in filepath_list:
        day_numb = (re.search(('Day_\d'), f)[0]).split('_')[-1]
        mouse_id = os.path.basename(f).split("_")[0]
        freq_arr = calc_frequency(f)
        freq_dict[mouse_id] = freq_arr
        if len(freq_dict.keys()) != len(filepath_list):
            pass
        else:
            df = pd.DataFrame.from_dict(freq_dict)
            # create dataframe from filepath and pair with event name
            group_behavior_df = (
                df
                .rolling(window=2, center=True).mean()
                .assign(
                    # mean=df.mean(axis=1),  # add mean,
                    time_sec=np.arange(-10, 21, 0.2),  # add time column
                    recording=recording,
                    day=(day_numb)
                )
                .dropna()
                .melt(id_vars=['day', 'time_sec', 'recording'], var_name='mouse', value_name='avg_frequency')
            )

    return recording, group_behavior_df


def collect_fp_data(path, sensor_regex, region_regex, event_regex):
    """function to collect photometry data and format for tidy format"""

    day_numb = (re.search(('Day_\d'), path)[0]).split('_')[-1]
    sensor = re.search(sensor_regex, os.path.basename(path))[0]
    region = re.search(region_regex, os.path.basename(path))[0]
    event = re.search(event_regex, os.path.basename(path))[0]

    df = pd.read_hdf(path)
    time_stamps = df['timestamps']
    df_clean = (
        df
        .drop(['timestamps'], axis=1)
        .rolling(window=500, center=True).mean()
        .assign(day=day_numb)
        .rename(columns=lambda c: c.split('-')[0])
        .drop(columns=['err'])
        .assign(time_stamps=time_stamps)
        .dropna()
        .melt(id_vars=['day', 'time_stamps'], var_name='mouse', value_name='z-score')
        .assign(sensor=sensor,
                region=region,
                event=event,
                recording='photometry')
    )
    return df_clean


def collect_AUC_data(path, sensor_regex, region_regex, event_regex):
    """function to collect AUC and peak data -> format for tidy format"""

    day_numb = (re.search(('Day_\d'), path)[0]).split('_')[-1]
    sensor = re.search(sensor_regex, os.path.basename(path))[0]
    region = re.search(region_regex, os.path.basename(path))[0]
    event = re.search(event_regex, os.path.basename(path))[0]

    df = pd.read_csv(path)
    AUC_clean = (
        df
        .assign(mouse=df['Unnamed: 0'].str.split('-', expand=True)[0],
                trial=df['Unnamed: 0'].str.split('_', expand=True)[1])
        .drop(columns='Unnamed: 0')
        .melt(id_vars=['mouse', 'trial'], var_name='category', value_name='value')
        .assign(sensor=sensor,
                region=region,
                event=event,
                day=day_numb,
                recording='photometry')
    )
    return AUC_clean


"""
GROUP AGGREGATION FUNCTIONS
"""


def compile_data(file_list):
    dataframes = [pd.read_feather(f) for f in file_list]
    compiled_df = pd.concat(dataframes)
    return compiled_df


"""
PLOTTING FUNCTIONS
"""


def fp_plot_line(df, sensor=None, region=None, event=None, sub_axes=None, alpha=None, color=None, linewidth=None):

    data = df[(df.sensor == sensor)
              & (df.region == region)
              & (df.event == event)
              ]

    sns.lineplot(data=data,
                 x='time_stamps',
                 y='z-score_mean',
                 color=color,
                 linewidth=linewidth,
                 ax=sub_axes
                 )
    sub_axes.fill_between(data=data, x=data['time_stamps'], y1=(data['z-score_mean'] + data['z-score_sem']),
                          y2=(data['z-score_mean'] - data['z-score_sem']), color=color, alpha=alpha, linewidth=0)


def behavior_plot_line(df, recording=None, sub_axes=None, alpha=None, color=None, linewidth=None, hue=None):

    data = df[(df.recording == recording)]

    sns.lineplot(data=data,
                 x='time_sec',
                 y='avg_frequency_mean',
                 color=color,
                 linewidth=linewidth,
                 hue=hue,
                 ax=sub_axes
                 )
    sub_axes.fill_between(data=data, x=data['time_sec'], y1=(data['avg_frequency_mean'] + data['avg_frequency_sem']),
                          y2=(data['avg_frequency_mean'] - data['avg_frequency_sem']), color=color, alpha=alpha, linewidth=0)


def draw_cue_box(ax, color, alpha):
    #  draw box on plot for cue
    y_lower = ax.get_ylim()[0]
    y_ags_sum = sum(np.abs(ax.get_ylim()))
    rect = patches.Rectangle(
        (0, y_lower), width=5, height=y_ags_sum, alpha=alpha, facecolor=color)
    ax.add_patch(rect)
    return rect
