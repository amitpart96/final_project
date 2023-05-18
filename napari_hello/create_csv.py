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
import scipy.stats as stats


# delete the segmentation images from the list of proteins
def delete_seg(patient, list_proteins):
    for f in os.listdir(patient):
        if os.path.isfile(os.path.join(patient, f)):
            if "Segmentation" in f:
                # print(f)
                list_proteins.remove(f)


# append data to the dataframe from region properties
def append_data(labels_max, props, str):
    row_data = []
    for index in range(0, labels_max):
        row_data.append(getattr(props[index], str))
    return row_data


# calculate the protein data for the cells
def protein_culc(list_proteins, patient, labels_max, props, df):
    for protein in list_proteins:
        protein_IMG = Image.open("{}\{}".format(patient, protein))
        print(f'protein calculate : {protein.split(".")[0]}')
        # convert img to np array
        protein_IMG = np.array(protein_IMG)
        col_pro = []
        for index in range(0, labels_max):
            list_of_indexes = getattr(props[index], 'coords')
            cell_size = getattr(props[index], 'area')
            x = protein_IMG[list_of_indexes].sum() / cell_size * 100
            col_pro.append(math.log(x + math.sqrt(1+math.pow(x, 2))))
        # calculate a z-score
        # z_score = stats.zscore(col_pro)
        df[protein.split(".")[0]] = col_pro
    # df.to_csv('csv.csv', index=False)
    return df


# create the labeled cell data matrix
def labeledcellData_matrix_create(shape):
    size1 = shape[0]
    size2 = shape[1]
    labeledcellData_matrix = np.ones((size1, size2))
    return labeledcellData_matrix


# update the labeled cell data matrix
def labeledcellData_matrix_update(labeledcellData_matrix, index, list_of_indexes):
    for x, y in list_of_indexes:
        labeledcellData_matrix[x, y] = index + 1
    return labeledcellData_matrix


# write the dataframe to a csv file.
def write_csv(file, df):
    if df is None:
        return "Error"
    df.to_csv(file.name, index=False)


# create the csv file
def create_csv(root_directory_path):
    global root_dir
    root = tk.Tk()
    root.withdraw()
    root_dir = root_directory_path
    print(f'root_dir: {root_dir}')
    # the user chooses the file name and the directory of the csv file
    file = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    if file is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    print(f'file name: {file.name}')
    return file


# This function finds the subfolders of the patients and returns a list of tuples with the path and name of each sub
# folder
def patient(path):
    # find the subfolders of the patients - each sujbfolder is one patient that contains his proteins and a segmantation
    list_subfolders_with_paths = [(f.path, f.name) for f in os.scandir(root_dir) if f.is_dir()]
    print(list_subfolders_with_paths)
    result = []

    # For each patient subfolder, create a dataframe with information about the cells in the segmentation image
    for patient in list_subfolders_with_paths:
        df = pd.DataFrame()
        print(f'patient: {patient[1]}')
        # Get a list of all the protein image files in the patient's subfolder
        list_proteins = ([f for f in os.listdir(patient[0]) if os.path.isfile(os.path.join(patient[0], f))])
        # filter the images of the proteins so that they will not contain the segmentation
        delete_seg(patient[0], list_proteins)

        # Load the segmentation image and use it to label the cells
        image = Image.open(patient[0] + '\SegmentationInterior.tiff')
        image = np.array(image)

        labels = measure.label(image, connectivity=2)
        props = measure.regionprops(labels)
        labels_max = labels.max()

        # Add columns to the dataframe with information about the cells
        df['SampleID'] = pd.Series([patient[1] for x in range(labels_max)])
        df.index = np.arange(1, len(df) + 1)
        df['cellLabelInImage'] = np.arange(2, len(df) + 2)
        col_sell_size = append_data(labels_max, props, 'area')
        df['cellSize'] = col_sell_size

        # Create a matrix to store information about each cell in the segmentation image
        labeledcellData_matrix = labeledcellData_matrix_create(image.shape)
        for index in range(0, labels_max):
            # Get the coordinates of each pixel in the cell and update the matrix with the cell's information
            list_of_indexes = getattr(props[index], 'coords')
            labeledcellData_matrix = labeledcellData_matrix_update(labeledcellData_matrix, index, list_of_indexes)

        # Save the matrix as an image file
        save_img(labeledcellData_matrix, "{}/p{}_{}".format(path,patient[1], "labeledcellData"))
        # Calculate information about each protein in the patient's subfolder and add it to the dataframe
        protein_culc(list_proteins, patient[0], labels_max, props, df)
        result.append(df)
        print(f"Done Calculate for patient: {patient[1]}")
    result = pd.concat(result)
    return result


# This function saves a matrix as an image file with the specified name
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


def main(root_directory_path):
    f = create_csv(root_directory_path)
    print(f)
    path = f.name
    path = path[:path.rfind('/')]
    print(f'file_path": {path}')
    data = patient(path)
    write_csv(f, data)
    print("created cellTable successfully")


#Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main(sys.argv)

