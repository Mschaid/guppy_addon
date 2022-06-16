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
class for plot timeseries data


this includes FP traces,
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
        for s in stores_list:
            if f str(os.listdir(self.output_path))):
                print(s)


mouse_1=Subject(r'R:\Mike\LHA_dopamine\LH_NAC_Headfix_FP\Photometry\Pav Training\2_color_pav\Sucrose_to_sucralose\Training\Day_9\512581-220531-110552\512581-220531-110552_output_1')
mouse_1.get_raw_data()




# %%
mouse_1=Subject(r'R:\Mike\LHA_dopamine\LH_NAC_Headfix_FP\Photometry\Pav Training\2_color_pav\Sucrose_to_sucralose\Training\Day_9\512581-220531-110552\512581-220531-110552_output_1')
dir=mouse_1.get_output_path()

files=[os.path(join(dir, file) for file in os.listdir(dir))]
print(files)

# def list_subdirs(directory):
#     """ given directory, returns list of all sub dirs as full path"""
#     return [os.path.join(directory, file) for file in os.listdir(directory)]

# list_subdirs(dir)


# %%
class Behavior:
    pass


class Viz:
    pass


mouse_1=Subject(r'R:\Mike\LHA_dopamine\LH_NAC_Headfix_FP\Photometry\Pav Training\2_color_pav\Sucrose_to_sucralose\Training\Day_9\512581-220531-110552\512581-220531-110552_output_1')


# %%
