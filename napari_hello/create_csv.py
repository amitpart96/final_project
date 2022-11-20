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




def create_csv():
    print("here0")
    all_rows = []
    col_names_csv = ['patient Number', 'cell index', 'cell size']
    all_rows.append(col_names_csv)

    root = tk.Tk()
    root.withdraw()

    # the current patient starts from index 1
    start_index_curr_patient = 1

    rootdir = filedialog.askdirectory()
    print("1")
    print(rootdir)
    print("2")

    # the user chooses the file name and the directory of the csv file
    file = filedialog.asksaveasfile(mode='w', defaultextension=".csv")
    if file is None:  # asksaveasfile return `None` if dialog closed with "cancel".
        return
    print(file.name)

    # find the subfolders of the patients - each subfolder is one patient that contains his proteins and a segmantation
    list_subfolders_with_paths = [f.path for f in os.scandir(rootdir) if f.is_dir()]
    for patient in list_subfolders_with_paths:
       # print(patient + '\GtSegmentationInterior.tiff')
       # print([f for f in os.listdir(patient) if os.path.isfile(os.path.join(patient, f))])
        list_proteins = ([f for f in os.listdir(patient) if os.path.isfile(os.path.join(patient, f))])
        # filter the images of the proteins so that they will not contain the segmentation
        for f in os.listdir(patient):
            if os.path.isfile(os.path.join(patient, f)) and f.endswith(".tiff"):
                if "Segmentation" in f:
                    print(f)
                    list_proteins.remove(f)
        print(list_proteins)
        # add to the column names of the csv with the proteins
        col_names_csv.extend(list_proteins)
        print(col_names_csv)
        all_rows.append(col_names_csv)

       # oldIMG = Image.open(patient + '\GtSegmentationInterior.tiff')

       # oldIMG = np.array(oldIMG)
        image = Image.open(patient + '\SegmentationInterior.tiff')
        image = np.array(image)

        # convert img to binary img
        #img = np.where(oldIMG == 0, 0, 1)

        # Binary image, post-process the binary mask and compute labels
        # threshold = filters.threshold_otsu(img)
        # mask = img > threshold
        # mask = morphology.remove_small_objects(mask, 100)
        # mask = morphology.remove_small_holes(mask, 100)
        # labels = measure.label(mask)
        # print(labels.max())

       # labels2 = measure.label(img)
       # print(labels2.max())
        labels = measure.label(image, connectivity=2)
       # print(labels.max())

        # print(labels3.num)

        # fig = px.imshow(img, binary_string=True)
        #  fig.update_traces(hoverinfo='skip')  # hover is only for label info
       # print(type(labels))

        props = measure.regionprops(labels)
       # properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

        for index in range(1, labels.max()):
            # label_i = props[index].label
            # contour = measure.find_contours(labels == label_i, 0.5)[0]
            # y, x = contour.T
            # hoverinfo = ''
            row_data = [patient] #To-Do: change to the folder name instead of the full path.
            row_data.append(index)
            row_data.append(getattr(props[index], 'area'))
            all_rows.append(row_data)

        # For each label, add a filled scatter trace for its contour,
        # and display the properties of the label in the hover of this trace.
        # for protein in proteins:




        for protein in list_proteins:
            # protein_IMG = Image.open("{}\{}.tiff".format(patient, protein))
            # print("{}\{}.tiff".format(patient, protein))
            protein_IMG = Image.open("{}\{}".format(patient, protein))
            print(protein)
            # print(protein_IMG)


            # convert img to np array
            protein_IMG = np.array(protein_IMG)
            for index in range(1, labels.max()):
                list_of_indexes = getattr(props[index], 'coords')

                #append to the current patient( we know the row of the current patient by adding the index of the cell to the start row of the current patient in the csv file
                all_rows[start_index_curr_patient + index - 1].append(protein_IMG[list_of_indexes].sum())

        #add the number of the cells of the current patient to the start index of "current patient" before moving to the next patient
        start_index_curr_patient += labels.max() - 1
        print(start_index_curr_patient)



    with open(file.name, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f)
        print("here1")
        # write the header
        writer.writerows(all_rows)
        print("here2")



def main():
    """
    oldIMG = Image.open('GtSegmentationInterior.tiff')
    #convert img to np array
    oldIMG = np.array(oldIMG)

    #convert img to binary img
    img = np.where(oldIMG ==0 , 0, 1)


    # Binary image, post-process the binary mask and compute labels
    threshold = filters.threshold_otsu(img)
    mask = img > threshold
    mask = morphology.remove_small_objects(mask, 100)
    mask = morphology.remove_small_holes(mask, 100)
    labels = measure.label(mask)

    fig = px.imshow(img, binary_string=True)
    fig.update_traces(hoverinfo='skip')  # hover is only for label info



    props = measure.regionprops(labels, img)
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']


    # For each label, add a filled scatter trace for its contour,
    # and display the properties of the label in the hover of this trace.
    for index in range(1, labels.max()):
        label_i = props[index].label
        contour = measure.find_contours(labels == label_i, 0.5)[0]
        y, x = contour.T
        hoverinfo = ''
        for prop_name in properties:
            hoverinfo += f'<b>{prop_name}: {getattr(props[index], prop_name):.2f}</b><br>'
        fig.add_trace(go.Scatter(
            x=x, y=y, name=label_i,
            mode='lines', fill='toself', showlegend=False,
            hovertemplate=hoverinfo, hoveron='points+fills'))


    plotly.io.show(fig)
    print("done show")
    """
    create_csv()
    cProfile.run('re.compile("foo|bar")')




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()




