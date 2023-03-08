from turtle import fd
from PIL.Image import Image
from imageio import imread
from napari_hello import ranking_model
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
import napari
from skimage import data


def find_anomaly(df, protein, patient):#todo: change name to find anomaly
    df = df.copy()
    print(f'starting patient number: {patient}')
    flag=True
    while flag:
        try:
            print("entered try")
            patient_labeled_cell_data = fd.askopenfilename()  # choose celldata of the patient
            print("1")
            cellLabel_image = Image.open(patient_labeled_cell_data)
            #cellLabel_image = Image.open('C:/Users/lidor/Downloads/Project2022/drive/TNBC_shareCellData/p1_labeledcellData.tiff')
            #from PIL import Image
            #todo: try to convert the image to list
            cellLabel_image = np.asarray(cellLabel_image)
            print("2")
            print("3")
            flag = False
            print("done try")
        except:
            print("incoreect path to celldata.tiff of the testing patient")

    print("before ranking model function call")
    DTR_scores, DTR_r2_scores, DTR_prediction = ranking_model.ranking_model(df, patient, protein)
    print("after ranking model function call")

    for protein, protein_prediction in DTR_prediction.items():  # DTR_prediction is a dictionary
        print(f'starting protein : {protein}')
        prediction_matrix = ranking_model.prediction_matrix_creation(protein_prediction, df, patient, cellLabel_image)
        # prediction matrix to image:
        pred_img=ranking_model.save_img(prediction_matrix, f'protein_prediction_{patient}_{protein}')
        real_protein_matrix = ranking_model.real_protein_matrix_creation(df, patient, protein, cellLabel_image)
        # real matrix to image:
        real_img=ranking_model.save_img(real_protein_matrix, f'real_protein_{patient}_{protein}')
        # std of the real matrix
        std_real = real_protein_matrix.std()
        print(std_real)
        # difference by std
        difference_matrix_std = create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real)
        difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
        # difference_matrix to image:
        diff_img = ranking_model.save_img(difference_matrix, f'difference_matrix_{patient}_{protein}')
        diff_img_std = ranking_model.save_img(difference_matrix_std, f'difference_matrix_std_{patient}_{protein}')
    print(f'finished patient number: {patient}')
    return real_img,pred_img,diff_img,diff_img_std

def create_difference_matrix_std(prediction_matrix,real_protein_matrix, std_real):
    difference_matrix_std_tmp = abs(np.subtract(prediction_matrix, real_protein_matrix))
    difference_matrix_std = (difference_matrix_std_tmp >= 2* std_real)
    print((f'std matrix Type: {type(difference_matrix_std)}'))
    return difference_matrix_std

def main(viewer,df,patient_number,protein):
    list_of_proteins_to_predict=[protein]
    real_img,pred_img,diff_img, diff_img_std=find_anomaly(df, list_of_proteins_to_predict, patient_number)
    # napari_image = imread(real_img)  # Reads an image from file
    # viewer.add_image(napari_image, name=real_img, colormap='magma')  # Adds the image to the viewer and give the image layer a name
    # napari_image = imread(pred_img)  # Reads an image from file
    # viewer.add_image(napari_image, name=pred_img, colormap='red')  # Adds the image to the viewer and give the image layer a name
    # napari_image = imread(diff_img)  # Reads an image from file
    # viewer.add_image(napari_image, name=diff_img, colormap='turbo')  # Adds the image to the viewer and give the image layer a name
    # napari_image = imread(diff_img_std)  # Reads an image from file
    # viewer.add_image(napari_image, name=diff_img)  # Adds the image to the viewer and give the image layer a nam


    # Add the real image as the first layer and set the colormap
    napari_image = imread(real_img)
    viewer.add_image(napari_image, name='Real Image', colormap='magma')

    # Add the predicted image as the second layer and set the colormap
    napari_image = imread(pred_img)
    viewer.add_image(napari_image, name='Predicted Image', colormap='red')

    # Add the difference image as the third layer and set the colormap
    napari_image = imread(diff_img)
    viewer.add_image(napari_image, name='Difference Image', colormap='turbo')

    # Add the difference image with std as the fourth layer and set the colormap
    napari_image = imread(diff_img_std)
    viewer.add_image(napari_image, name='Difference Image with Std', colormap='gray')

    # Add the labels image as the fifth layer and set the colormap
    labels_image = imread(patient_labeled_cell_data)
    viewer.add_labels(labels_image, name='Labels Image', colormap='viridis')

if __name__ == "__main__":
    main()
