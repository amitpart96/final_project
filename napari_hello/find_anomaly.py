from turtle import fd
from PIL.Image import Image
from imageio import imread
from napari_hello import ranking_model
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image


def find_anomaly(df, protein, patient, model_name,proteins_list,patient_cellLabel_image, slider_float):  # todo: change name to find anomaly
    """
    Find anomalies in protein predictions and generate visualizations.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing protein data.
    - protein (str): The name of the protein to analyze.
    - patient (str): The identifier of the patient.
    - model_name (str): The name of the predicion model.
    - proteins_list (list): A list of proteins to consider in the analysis.
    - patient_cellLabel_image (np.ndarray): The cell label image for the patient.

    Returns:
    Tuple: A tuple containing the following:
    - real_img (str): The filename of the saved image showing the real protein distribution.
    - pred_img (str): The filename of the saved image showing the predicted protein distribution.
    - diff_img (str): The filename of the saved image showing the difference between the real and predicted protein distributions.
    - diff_img_std (str): The filename of the saved image showing the difference between the real and predicted protein distributions normalized by standard deviation.
    - prediction_matrix (np.ndarray): The matrix representing the predicted protein distribution.
    - real_protein_matrix (np.ndarray): The matrix representing the real protein distribution.
    - std_real (float): The standard deviation of the real protein distribution.
    - file_name_std (str): The filename of the saved image showing the difference matrix normalized by standard deviation.
    """
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
        # difference by std
        difference_matrix_std = create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real, slider_float)
        difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
        # difference_matrix to image:
        diff_img = ranking_model.save_img(difference_matrix, f'difference_matrix_{patient}_{protein}')
        diff_img_std = ranking_model.save_img(difference_matrix_std, f'difference_matrix_std_{patient}_{protein}')
        file_name_std = f'difference_matrix_std_{patient}_{protein}'
    return real_img, pred_img, diff_img, diff_img_std, prediction_matrix, real_protein_matrix, std_real, file_name_std  # ??


def create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real, slider_float):
    """
    Create a difference matrix based on the standard deviation.

    Parameters:
    - prediction_matrix (np.ndarray): The matrix representing the predicted protein distribution.
    - real_protein_matrix (np.ndarray): The matrix representing the real protein distribution.
    - std_real (float): The standard deviation of the real protein distribution.

    Returns:
    np.ndarray: The difference matrix where values greater than or equal to 2 times the standard deviation are True,
                and values less than 2 times the standard deviation are False.
    """
    difference_matrix_std_tmp = abs(np.subtract(prediction_matrix, real_protein_matrix))
    difference_matrix_std = (difference_matrix_std_tmp >= slider_float * std_real)
    return difference_matrix_std


def update_difference_matrix_std(viewer, prediction_matrix, real_protein_matrix, std_real, slider_float, file_name_std,
                                 layer_std):
    """
    Update the difference matrix std based on the standard deviation and display it in the viewer.

    Parameters:
    - viewer (napari.viewer.Viewer): The Napari viewer object.
    - prediction_matrix (np.ndarray): The matrix representing the predicted protein distribution.
    - real_protein_matrix (np.ndarray): The matrix representing the real protein distribution.
    - std_real (float): The standard deviation of the real protein distribution.
    - slider_float (float): The slider value representing the threshold multiplier.
    - file_name_std (str): The file name of the saved difference matrix image.
    - layer_std (napari.layers.Image): The Napari image layer for the difference matrix.

    Returns:
    None
    """
    difference_matrix_std_tmp = abs(np.subtract(prediction_matrix, real_protein_matrix))
    difference_matrix_std = (difference_matrix_std_tmp >= slider_float * std_real)
    diff_img_std = ranking_model.save_img(difference_matrix_std, file_name_std)
    napari_image = imread(diff_img_std)  # Reads an image from file
    # viewer.add_image(napari_image, name=diff_img_std)  # Adds the image to the viewer and give the image layer a name
    # viewer.layers[file_name_std].data = diff_img_std #'(napari_image, name=file_name_std)  # Adds the image to the viewer and give the image layer a name
    layer_std.data = napari_image


def main(viewer, df, patient_number, protein, model_name,proteins_list,patient_cellLabel_image, slider_float):
    list_of_proteins_to_predict = [protein]
    real_img, pred_img, diff_img, diff_img_std, prediction_matrix, real_protein_matrix, std_real, file_name_std = find_anomaly(
        df, list_of_proteins_to_predict, patient_number, model_name,proteins_list,patient_cellLabel_image, slider_float)  # ??
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
