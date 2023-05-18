import os
import sys
from turtle import fd
from PIL.Image import Image
from imageio import imread
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image

import ranking_model


def find_anomaly(df, protein, patient, model_name,proteins_list,patient_cellLabel_image,save_path,std):
    df = df.copy()
    print(f'starting patient number: {patient}')
    print("before ranking model function call")
    scores, r2_scores, prediction = ranking_model.ranking_model(df, patient, protein,proteins_list, model_name)
    print("after ranking model function call")

    for protein, protein_prediction in prediction.items():  # DTR_prediction is a dictionary
        print(f'protein_prediction: {protein_prediction}')
        print(f'starting protein : {protein}')
        prediction_matrix = ranking_model.prediction_matrix_creation(protein_prediction, df, patient, patient_cellLabel_image)
        # prediction matrix to image:
        pred_img = ranking_model.save_img(prediction_matrix,"{}/protein_prediction{}_{}".format(save_path, patient, protein))
        real_protein_matrix = ranking_model.real_protein_matrix_creation(df, patient, protein, patient_cellLabel_image)
        # real matrix to image:
        real_img = ranking_model.save_img(real_protein_matrix,"{}/real_protein{}_{}".format(save_path, patient, protein))
        # std of the real matrix
        std_real = real_protein_matrix.std()
        print(std_real)
        # difference by std
        difference_matrix_std = create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real, std)
        difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
        # difference_matrix to image:
        diff_img = ranking_model.save_img(difference_matrix,"{}/difference_matrix_{}_{}".format(save_path, patient, protein))
        diff_img_std = ranking_model.save_img(difference_matrix_std,"{}/difference_matrix_std_{}_{}".format(save_path, patient, protein))
        file_name_std = f'difference_matrix_std_{patient}_{protein}'
    print(f'finished patient number: {patient}')
    return real_img, pred_img, diff_img, diff_img_std, prediction_matrix, real_protein_matrix, std_real, file_name_std


def create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real, std):
    difference_matrix_std_tmp = abs(np.subtract(prediction_matrix, real_protein_matrix))
    difference_matrix_std = (difference_matrix_std_tmp >= std * std_real)
    print((f'std matrix Type: {type(difference_matrix_std)}'))
    return difference_matrix_std


def main(both ,df, patient_number, protein, model_name, save_path,std):

    if both != 'True':
        df = pd.read_csv(df)

    path_imag = "{}\p{}_{}".format(save_path, patient_number, "labeledcellData.tiff")
    cellLabel_image = Image.open(path_imag)
    array_img = np.asarray(cellLabel_image)
    print(array_img)
    column_names = df.columns.tolist()
    # Columns to remove
    columns_to_remove = ['SampleID', 'cellLabelInImage', 'cellSize']  # List of column names to remove
    proteins_list = [col for col in column_names if col not in columns_to_remove]
    print(df)
    list_of_proteins_to_predict = [protein]
    find_anomaly(df, list_of_proteins_to_predict, patient_number, model_name,proteins_list,array_img,save_path,std)

