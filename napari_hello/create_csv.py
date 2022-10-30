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



def create_csv():
    print("here0")
    all_rows = []
    col_names_csv = ['patient Number', 'cell index', 'cell size', 'Au','B7H3', 'Background', 'Beta catenin', 'C', 'Ca', 'CD3', 'CD4', 'CD8', 'CD11b', 'CD11c', 'CD16', 'CD20' , 'CD31', 'CD45', 'CD45RO', 'CD56', 'CD63', 'CD68', 'CD68', 'CD163', 'CD209', 'CD209', 'CellTypes', 'CSF-1R', 'dsDNA','EGFR','Fe','FoxP3','H3K9ac','H3K27me3','HLA_class_1','HLA-DR','IDO','Keratin6','Keratin17','Ki67','Lag3','MPO','Na','OX40','p','p53','Pan-Keratin','PD1','PD-L1','phospho-S6','Si','SMA','Ta','Vimentin' ]
    #list of all the proteins
    #proteins = ['Au', 'B7H3','Background', 'Beta catenin', 'C', 'Ca', 'CD3', 'CD4', 'CD8', 'CD11b', 'CD11c', 'CD16', 'CD20' , 'CD31', 'CD45', 'CD45RO', 'CD56', 'CD63', 'CD68', 'CD68', 'CD163', 'CD209', 'CD209', 'CellTypes', 'CSF-1R','dsDNA','EGFR','Fe','FoxP3','H3K9ac','H3K27me3','HLA_class_1','HLA-DR','IDO','Keratin6','Keratin17','Ki67','Lag3','MPO','Na','OX40','p','p53','Pan-Keratin','PD1','PD-L1','phospho-S6','Si','SMA','Ta','Vimentin']
    proteins = ['Au', 'B7H3']
    all_rows.append(col_names_csv)

    root = tk.Tk()
    root.withdraw()



    rootdir  = filedialog.askdirectory()
    print("1")
    print(rootdir)
    print("2")
    list_subfolders_with_paths = [f.path for f in os.scandir(rootdir) if f.is_dir()]
    for patient in list_subfolders_with_paths:
        print(patient+'\GtSegmentationInterior.tiff')

        # convert img to np array
        oldIMG = np.array(Image.open(patient+'\GtSegmentationInterior.tiff'))

        # convert img to binary img
        img = np.where(oldIMG == 0, 0, 1)

        # Binary image, post-process the binary mask and compute labels
        threshold = filters.threshold_otsu(img)
        mask = img > threshold
        mask = morphology.remove_small_objects(mask, 100)
        mask = morphology.remove_small_holes(mask, 100)
        labels = measure.label(mask)

        # fig = px.imshow(img, binary_string=True)
        # fig.update_traces(hoverinfo='skip')  # hover is only for label info

        props = measure.regionprops(labels, img)
        # properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']

        # For each label, add a filled scatter trace for its contour,
        # and display the properties of the label in the hover of this trace.
        for protein in proteins:
            protein_IMG = Image.open("{}\{}.tiff".format(patient, protein))
            print("{}\{}.tiff".format(patient, protein))
            #print(protein_IMG)
            # convert img to np array

            protein_IMG = np.array(protein_IMG)
            for index in range(1, labels.max()):
                # print(index)
                label_i = props[index].label
                contour = measure.find_contours(labels == label_i, 0.5)[0]
                y, x = contour.T
                hoverinfo = ''
                row_data = [patient]
                row_data.append(index)
                row_data.append(getattr(props[index], 'area'))
                list_of_indexes = getattr(props[index], 'coords')


                row_data.append(protein_IMG[list_of_indexes].sum())
                all_rows.append(row_data)
        # for prop_name in properties:
        #
        #     hoverinfo += f'<b>{prop_name}: {getattr(props[index], prop_name):.2f}</b><br>'

        # fig.add_trace(go.Scatter(
        #     x=x, y=y, name=label_i,
        #     mode='lines', fill='toself', showlegend=False,
        #     hovertemplate=hoverinfo, hoveron='points+fills'))
    # plotly.io.show(fig)
    # print("done show")

    with open('resultNew2.csv', 'w', encoding='UTF8', newline='') as f:
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



#
# # Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     main()




