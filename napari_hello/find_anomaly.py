from turtle import fd
from PIL.Image import Image
from imageio import imread
from napari_hello import ranking_model
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image


def find_anomaly(df, protein, patient, model_name,proteins_list,patient_cellLabel_image):  # todo: change name to find anomaly
    df = df.copy()
    scores, r2_scores, prediction = ranking_model.ranking_model(df, patient, protein,proteins_list, model_name)

    for protein, protein_prediction in prediction.items():  # DTR_prediction is a dictionary
        prediction_matrix = ranking_model.prediction_matrix_creation(protein_prediction, df, patient, patient_cellLabel_image)
        # prediction matrix to image:
        pred_img = ranking_model.save_img(prediction_matrix, f'protein_prediction_{patient}_{protein}')
        real_protein_matrix = ranking_model.real_protein_matrix_creation(df, patient, protein, patient_cellLabel_image)
        # real matrix to image:
        real_img = ranking_model.save_img(real_protein_matrix, f'real_protein_{patient}_{protein}')
        # std of the real matrix
        std_real = real_protein_matrix.std()
        print(std_real)
        # difference by std
        difference_matrix_std = create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real)
        difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
        # difference_matrix to image:
        diff_img = ranking_model.save_img(difference_matrix, f'difference_matrix_{patient}_{protein}')
        diff_img_std = ranking_model.save_img(difference_matrix_std, f'difference_matrix_std_{patient}_{protein}')
        file_name_std = f'difference_matrix_std_{patient}_{protein}'
    return real_img, pred_img, diff_img, diff_img_std, prediction_matrix, real_protein_matrix, std_real, file_name_std  # ??


def create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real):
    difference_matrix_std_tmp = abs(np.subtract(prediction_matrix, real_protein_matrix))
    difference_matrix_std = (difference_matrix_std_tmp >= 2 * std_real)
    return difference_matrix_std


def update_difference_matrix_std(viewer, prediction_matrix, real_protein_matrix, std_real, slider_float, file_name_std,
                                 layer_std):
    difference_matrix_std_tmp = abs(np.subtract(prediction_matrix, real_protein_matrix))
    difference_matrix_std = (difference_matrix_std_tmp >= slider_float * std_real)
    diff_img_std = ranking_model.save_img(difference_matrix_std, file_name_std)
    napari_image = imread(diff_img_std)  # Reads an image from file
    # viewer.add_image(napari_image, name=diff_img_std)  # Adds the image to the viewer and give the image layer a name
    # viewer.layers[file_name_std].data = diff_img_std #'(napari_image, name=file_name_std)  # Adds the image to the viewer and give the image layer a name
    layer_std.data = napari_image


def main(viewer, df, patient_number, protein, model_name,proteins_list,patient_cellLabel_image):
    list_of_proteins_to_predict = [protein]
    real_img, pred_img, diff_img, diff_img_std, prediction_matrix, real_protein_matrix, std_real, file_name_std = find_anomaly(
        df, list_of_proteins_to_predict, patient_number, model_name,proteins_list,patient_cellLabel_image)  # ??
    napari_image = imread(real_img)  # Reads an image from file
    viewer.add_image(napari_image, name=real_img)  # Adds the image to the viewer and give the image layer a name
    napari_image = imread(pred_img)  # Reads an image from file
    viewer.add_image(napari_image, name=pred_img)  # Adds the image to the viewer and give the image layer a name
    napari_image = imread(diff_img)  # Reads an image from file
    viewer.add_image(napari_image, name=diff_img)  # Adds the image to the viewer and give the image layer a name
    napari_image = imread(diff_img_std)  # Reads an image from file
    layer_std = viewer.add_image(napari_image,
                                 name=diff_img_std)  # Adds the image to the viewer and give the image layer a name
    return prediction_matrix, real_protein_matrix, std_real, file_name_std, layer_std


if __name__ == "__main__":
    main()
