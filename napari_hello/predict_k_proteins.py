from collections import Counter
from random import randint
import random
import pandas as pd
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from skimage.io import imread
import time
import xgboost as xgb
from napari_hello import ranking_model


def predict_k_proteins(viewer, df, patient_number, list_of_proteins_to_predict, proteins_list, model_name,cellLabel_image):
    """
    Predict the values of multiple proteins for a given patient and visualize the predictions in the viewer.

    Parameters:
    - viewer (napari.viewer.Viewer): The Napari viewer object.
    - df (pd.DataFrame): The DataFrame containing the protein data.
    - patient_number (int): The identifier of the patient for prediction.
    - list_of_proteins_to_predict (list): A list of protein names to predict.
    - proteins_list (list): A list of all available protein names.
    - model_name (str): The name of the prediction model to use (either 'XGBoost' or 'DecisionTree').
    - cellLabel_image (np.ndarray): The cell label image corresponding to the patient.

    Returns:
    - DTR_cor_score (float): The correlation score of the predictions.
    - DTR_r2_score (float): The R2 score of the predictions.
    - prediction (np.ndarray): The predicted protein values.
    """
    df = df.copy()
    df_train = df.loc[df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
    df_test = df.loc[df['SampleID'] == patient_number]  # takes only patient patient_number for test

    # predict one protein , we will put it inside Y_train:
    if (len(list_of_proteins_to_predict) == 1):
        y_train, y_test = df_train[list_of_proteins_to_predict[0]], df_test[list_of_proteins_to_predict[0]]
    else:
        y_train, y_test = df_train[list_of_proteins_to_predict], df_test[list_of_proteins_to_predict]

    # we will put all the rest proteins inside X_train:
    pl_copy = proteins_list.copy()
    for protein in list_of_proteins_to_predict:
        pl_copy.remove(protein)

    X_train, X_test = df_train[pl_copy], df_test[pl_copy]

    if (model_name == 'XGBoost'):
        DTR_cor_score, DTR_r2_score, prediction = ranking_model.model_XGBoostRegressor(X_train, y_train, X_test, y_test)
    else:
        DTR_cor_score, DTR_r2_score, prediction = ranking_model.model_DecisionTreeRegressor(X_train, y_train, X_test, y_test)


    cor_plt_name,r2_plt_name= check_prediction_score(df, df_test, patient_number, prediction, list_of_proteins_to_predict, proteins_list,
                           model_name)

    prediction_df = pd.DataFrame(prediction, columns=list_of_proteins_to_predict)

    for protein_name, values in prediction_df.iteritems():
        protein_prediction_matrix = ranking_model.prediction_matrix_creation(prediction_df[protein_name], df, patient_number,
                                                               cellLabel_image)

        img_name = f'protein_prediction_{patient_number}_{protein_name}'
        img = ranking_model.save_img(protein_prediction_matrix, img_name)
        protein_prediction_image = imread(img)
        viewer.add_image(protein_prediction_image, name=img_name)  # Adds the image to the viewer and give the image layer a name

    viewer.add_image(cor_plt_name, name="Correlation comparison")
    viewer.add_image(r2_plt_name, name="r2 comparison")
    return DTR_cor_score, DTR_r2_score, prediction

def check_prediction_score(df, df_test, patient_number, multy_prediction, list_of_proteins_to_predict, proteins_list, model_name):
    """
    Compare the prediction scores for multiple proteins and generate visualizations.

    Parameters:
    - df (pd.DataFrame): The DataFrame containing the protein data.
    - df_test (pd.DataFrame): The DataFrame containing the test data for the patient.
    - patient_number (int): The identifier of the patient for prediction.
    - multy_prediction (np.ndarray): The predicted protein values for multiple proteins.
    - list_of_proteins_to_predict (list): A list of protein names to predict.
    - proteins_list (list): A list of all available protein names.
    - model_name (str): The name of the prediction model used.

    Returns:
    - cor_plt_name (str): The filename of the correlation comparison plot.
    - r2_plt_name (str): The filename of the R2 score comparison plot.
    """
    cor_scores_counter = Counter()
    r2_scores_counter = Counter()
    cor_scores, r2_scores, single_predictions = ranking_model.ranking_model(df, patient_number, list_of_proteins_to_predict,
                                                              proteins_list, model_name)

    # Reshape multy_prediction to 2D array if it has only one column
    if multy_prediction.ndim == 1:
        multy_prediction = multy_prediction.reshape(-1, 1)

    # Loop over each column
    for col_idx, protein_name in enumerate(list_of_proteins_to_predict):
        real_protein_values = df_test[protein_name]

        # check scores for multi protein prediction
        protein_multi_prediction = multy_prediction[:, col_idx]
        multi_prediction_cor_score = ranking_model.calculate_correlation(real_protein_values, protein_multi_prediction)
        multi_prediction_r2_score = ranking_model.calculate_r2_score(real_protein_values, protein_multi_prediction)
        # check scores for one protein prediction
        single_prediction_cor_score = cor_scores[protein_name]
        single_prediction_r2_score = r2_scores[protein_name]
        # save in counter
        cor_scores_counter[protein_name] = [multi_prediction_cor_score[0, 1], single_prediction_cor_score]
        r2_scores_counter[protein_name] = [multi_prediction_r2_score, single_prediction_r2_score]
    cor_plt_name= plot_comparison_prediction_scores(cor_scores_counter, "Correlation")
    r2_plt_name=plot_comparison_prediction_scores(r2_scores_counter, "r2")
    return cor_plt_name,r2_plt_name

def plot_comparison_prediction_scores(cor_scores_counter,score_name):
    """
    Generate a comparison plot for prediction scores.

    Parameters:
    - cor_scores_counter (Counter): A counter object containing the prediction scores.
    - score_name (str): The name of the prediction score (e.g., "Correlation", "R2").

    Returns:
    - napari_image (np.ndarray): The generated comparison plot image.
    """

    # Extract protein names and scores
    protein_names = list(cor_scores_counter.keys())
    scores1 = [cor_scores_counter[protein_name][0] for protein_name in protein_names]
    scores2 = [cor_scores_counter[protein_name][1] for protein_name in protein_names]

    # Convert protein_names to numeric indices
    x_indices = np.arange(len(protein_names))
    plt.clf()

    # Set up bar plot
    plt.bar(x_indices, scores1, width=0.35, label='multy prediction')
    plt.bar(x_indices + 0.35, scores2, width=0.35, label='single prediction')

    # Set x-ticks to protein names
    plt.xticks(x_indices + 0.35 / 2, protein_names)

    # Add labels and legend
    plt.xlabel('Protein')
    plt.ylabel(score_name + ' Scores')
    plt.title('Comparison of Prediction Scores')
    plt.legend()

    # the plot
    # plt.show()
    plt.savefig(score_name+"_comparison", dpi=170)
    napari_image= imread(score_name+"_comparison"+".png")  # Reads an image from file
    return napari_image
if __name__ == "__main__":
    predict_k_proteins(viewer, df, patient_number, list_of_proteins_to_predict, proteins_list, model_name, cellLabel_image)


