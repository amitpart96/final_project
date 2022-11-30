import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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
def delete_seg(patient, list_proteins):
    for f in os.listdir(patient):
        if os.path.isfile(os.path.join(patient, f)) and f.endswith(".tiff"):
            if "Segmentation" in f:
                print(f)
                list_proteins.remove(f)

def append_data(labels_max, props,str):
    row_data=[]
    for index in range(0, labels_max):
        # To-Do: change to the folder name instead of the full path.
        row_data.append(getattr(props[index], str))
    return row_data


def protein_culc(list_proteins, patient, labels_max, props,df):
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

def write_csv(file, df):
    if df is None:
        return ("Erorr")
    df.to_csv(file.name)



def create_csv():
    global root_dir
    print("here0")
    root = tk.Tk()
    root.withdraw()
    root_dir = filedialog.askdirectory()
    print("1")
    print(root_dir)
    print("2")
    # the user chooses the file name and the directory of the csv file
    file = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    if file is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    print(file.name)
    return file



def patient():

    # find the subfolders of the patients - each subfolder is one patient that contains his proteins and a segmantation
    list_subfolders_with_paths = [f.path for f in os.scandir(root_dir) if f.is_dir()]
    print(list_subfolders_with_paths)
    result = []

    for patient in list_subfolders_with_paths:
        df = pd.DataFrame()
        print(patient)

        list_proteins = ([f for f in os.listdir(patient) if os.path.isfile(os.path.join(patient, f))])
        # filter the images of the proteins so that they will not contain the segmentation
        delete_seg(patient, list_proteins)

        image = Image.open(patient + '\SegmentationInterior.tiff')
        image = np.array(image)
        labels = measure.label(image, connectivity=2)
        props = measure.regionprops(labels)
        labels_max = labels.max()

        df['patient Number'] = pd.Series([patient for x in range(labels_max)])
        df.index = np.arange(1, len(df)+1)
        df['cell index'] = np.arange(1,len(df)+1)
        col_sell_size = append_data(labels_max, props,'area')
        df['cell_size'] = col_sell_size


        protein_culc(list_proteins, patient, labels_max, props,df)
        result.append(df)
    result = pd.concat(result)
    return result



def main():
    # get the start time
    st = time.time()

    f = create_csv()
    print(f)
    data = patient()
    write_csv(f,data)

    # get the end time
    et = time.time()

    # get the execution time
    elapsed_time = et - st
    print('Execution time:', elapsed_time, 'seconds')
    final_res = elapsed_time / 60
    print('Execution time:', final_res, 'minutes')

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()




