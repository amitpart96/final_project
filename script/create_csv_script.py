# This is a sample Python script.
import sys
import math
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
import re


# def append_data(labels_max, props, str):
#     row_data = []
#     for index in range(0, labels_max):
#         row_data.append(getattr(props[index], str))
#     return row_data
#
#
# def protein_culc(list_proteins, patient, labels_max, props, df):
#     for protein in list_proteins:
#         protein_IMG = Image.open("{}\{}".format(patient, protein))
#         print(protein)
#         # convert img to np array
#         protein_IMG = np.array(protein_IMG)
#         col_pro = []
#         for index in range(0, labels_max):
#             list_of_indexes = getattr(props[index], 'coords')
#             col_pro.append(protein_IMG[list_of_indexes].sum())
#
#         df[protein] = col_pro
#         df.to_csv('csv.csv')
#     return df

# def create_csv(argv):
#     file = argv[1]
#     save_file_name = argv[2]
#     list_subfolders_with_paths = [f.path for f in os.scandir(file) if f.is_dir()]
#     result = []
#     for patient in list_subfolders_with_paths:
#         list_proteins = ([f for f in os.listdir(patient) if os.path.isfile(os.path.join(patient, f))])
#         for f in os.listdir(patient):
#             if os.path.isfile(os.path.join(patient, f)) and f.endswith(".tiff"):
#                 if "Segmentation" in f:
#                     list_proteins.remove(f)
#         image = Image.open(patient + '\SegmentationInterior.tiff')
#         image = np.array(image)
#         labels = measure.label(image, connectivity=2)
#         props = measure.regionprops(labels)
#         labels_max = labels.max()
#         df = pd.DataFrame()
#         path = patient.split('\\')
#         patient_number = path[len(path) - 1]
#         print(patient_number)
#         df['patient Number'] = pd.Series([patient_number for x in range(labels_max)])
#         df.index = np.arange(1, len(df) + 1)
#         df['cell index'] = np.arange(1, len(df) + 1)
#         col_sell_size = append_data(labels_max, props, 'area')
#         df['cell_size'] = col_sell_size
#         protein_culc(list_proteins, patient, labels_max, props, df)
#         result.append(df)
#     result = pd.concat(result)
#     print(save_file_name)
#     result.to_csv(save_file_name)
#
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
        print(protein.split(".")[0])
        # convert img to np array
        protein_IMG = np.array(protein_IMG)
        col_pro = []
        for index in range(0, labels_max):
            list_of_indexes = getattr(props[index], 'coords')
            cell_size = getattr(props[index], 'area')
            x = protein_IMG[list_of_indexes].sum()/cell_size * 100
            col_pro.append(math.log(x + math.sqrt(1+math.pow(x,2))))

        df[protein.split(".")[0]] = col_pro
        df.to_csv('csv.csv', index = False)
    return df

def labeledcellData_matrix_create(patient):
    labeledcellData_matrix = np.zeros((2048, 2048))
    return labeledcellData_matrix


def labeledcellData_matrix_update(labeledcellData_matrix, index, list_of_indexes):
    for x,y in list_of_indexes:
        labeledcellData_matrix[x,y] = index+1

    return labeledcellData_matrix


def patient(argv):
    file = argv[1]
    save_file_name = argv[2]
    # find the subfolders of the patients - each sujbfolder is one patient that contains his proteins and a segmantation
    list_subfolders_with_paths = [(f.path,f.name) for f in os.scandir(file) if f.is_dir()]
    print(list_subfolders_with_paths)
    result = []

    for patient in list_subfolders_with_paths:
        df = pd.DataFrame()
        print(patient)

        list_proteins = ([f for f in os.listdir(patient[0]) if os.path.isfile(os.path.join(patient[0], f))])
        # filter the images of the proteins so that they will not contain the segmentation
        delete_seg(patient[0], list_proteins)

        image = Image.open(patient[0] + '\SegmentationInterior.tiff')
        image = np.array(image)
        labels = measure.label(image, connectivity=2)
        props = measure.regionprops(labels)
        labels_max = labels.max()

        df['SampleID'] = pd.Series([patient[1] for x in range(labels_max)])
        df.index = np.arange(1, len(df)+1)
        df['cellLabelInImage'] = np.arange(1,len(df)+1)
        col_sell_size = append_data(labels_max, props,'area')


        df['cellSize'] = col_sell_size

        labeledcellData_matrix = labeledcellData_matrix_create(patient[0])
        for index in range(0, labels_max):
            list_of_indexes = getattr(props[index], 'coords')
            labeledcellData_matrix = labeledcellData_matrix_update(labeledcellData_matrix, index, list_of_indexes)

        save_img(labeledcellData_matrix,"{}_{}".format(patient[0], "labeledcellData"))
        protein_culc(list_proteins, patient[0], labels_max, props,df)
        result.append(df)
    result = pd.concat(result)
    result.to_csv(save_file_name)
    print("Done")
    return result

def save_img(matrix, file_name):
    matrix = matrix.astype(np.uint16)
    print(np.amax(matrix))
    new_im = Image.fromarray(matrix)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    print(image_filename)
    # save image using extension
    new_im.save(image_filename)
    return image_filename


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    # create_csv(sys.argv)
    patient(sys.argv)


