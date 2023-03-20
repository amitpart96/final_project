import math
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
import csv
import plotly
import plotly.express as px
import plotly.graph_objects as go
from skimage import data, filters, measure, morphology
import tkinter as tk
from tkinter import filedialog
import os
import cProfile
import re
import time
import pandas as pd
import sys


def append_data(labels_max, props, str):
    row_data = []
    for index in range(0, labels_max):
        row_data.append(getattr(props[index], str))
    return row_data


def protein_culc(list_proteins, patient, labels_max, props, df):
    for protein in list_proteins:
        protein_IMG = Image.open("{}\{}".format(patient, protein))
        print(protein)
        # convert img to np array
        protein_IMG = np.array(protein_IMG)
        col_pro = []
        for index in range(0, labels_max):
            list_of_indexes = getattr(props[index], 'coords')
            col_pro.append(protein_IMG[list_of_indexes].sum())

        df[protein] = col_pro
        df.to_csv('csv.csv')
    return df


def create_csv(argv):
    file = argv[1]
    save_file_name = argv[2]
    # inputFileName = input("Enter name of input file: ")
    list_subfolders_with_paths = [f.path for f in os.scandir(file) if f.is_dir()]
    result = []
    for patient in list_subfolders_with_paths:
        list_proteins = ([f for f in os.listdir(patient) if os.path.isfile(os.path.join(patient, f))])
        for f in os.listdir(patient):
            if os.path.isfile(os.path.join(patient, f)) and f.endswith(".tiff"):
                if "Segmentation" in f:
                    list_proteins.remove(f)
        image = Image.open(patient + '\SegmentationInterior.tiff')
        image = np.array(image)
        labels = measure.label(image, connectivity=2)
        props = measure.regionprops(labels)
        labels_max = labels.max()
        df = pd.DataFrame()
        path = patient.split('\\')
        patient_number = path[len(path) - 1]
        print(patient_number)
        df['patient Number'] = pd.Series([patient_number for x in range(labels_max)])
        df.index = np.arange(1, len(df) + 1)
        df['cell index'] = np.arange(1, len(df) + 1)
        col_sell_size = append_data(labels_max, props, 'area')
        df['cell_size'] = col_sell_size
        protein_culc(list_proteins, patient, labels_max, props, df)
        result.append(df)
    result = pd.concat(result)
    result.to_csv(save_file_name)

    
# Press the green button in the gutter to run the script.
if __name__ == '__main__':
     create_csv(sys.argv)
    # test()
