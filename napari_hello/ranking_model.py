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


def save_img(matrix, file_name):
    """
    Saves the matrix as a TIFF image file with the specified file name.

    Parameters:
    - matrix: The matrix as a NumPy array.
    - file_name: The desired file name for the saved image.

    Returns:
    - image_filename: The name of the saved image file.
    """
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    image_filename = f'{file_name}.tiff'
    new_im.save(image_filename)
    return image_filename


def ranking_model(df, patient_number, list_of_proteins_to_predict, proteins_list, model_name):
    """
    Performs ranking and prediction for a list of proteins using a specified model.

    Parameters:
    - df: The input DataFrame containing protein data.
    - patient_number: The patient number for whom predictions are made.
    - list_of_proteins_to_predict: A list of proteins to predict.
    - proteins_list: The complete list of proteins in the DataFrame.
    - model_name: The name of the model to use for prediction ('XGBoost' or Decision-Tree).

    Returns:
    - DTR_cor_scores: A Counter object containing correlation scores for each predicted protein.
    - DTR_r2_scores: A Counter object containing R-squared scores for each predicted protein.
    - DTR_predictions: A Counter object containing predictions for each protein.
    """
    df = df.copy()
    DTR_cor_scores, DTR_r2_scores, DTR_predictions = Counter(), Counter(), Counter()

    df_train = df.loc[
        df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
    df_test = df.loc[df['SampleID'] == patient_number]  # takes only patient patient_number for test

    for protein in list_of_proteins_to_predict:
        # predict one protein , we will put it inside Y_train:
        y_train, y_test = df_train[protein], df_test[protein]
        # we will put all the rest proteins inside X_train:
        pl_copy = proteins_list.copy()
        pl_copy.remove(protein)
        X_train, X_test = df_train[pl_copy], df_test[pl_copy]

        if (model_name == 'XGBoost'):
            DTR_cor_score, DTR_r2_score, DTR_prediction = model_XGBoostRegressor(X_train, y_train, X_test, y_test)
        else:
            DTR_cor_score, DTR_r2_score, DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test, y_test)

        DTR_cor_scores[protein] = float(DTR_cor_score[0, 1])
        DTR_r2_scores[protein] = DTR_r2_score
        DTR_predictions[protein] = DTR_prediction

    return DTR_cor_scores, DTR_r2_scores, DTR_predictions


def ranking_model_avg(df, list_of_proteins_to_predict, proteins_list, model_name, amount_of_patients):
    """
    Performs ranking and prediction for a list of proteins using a specified model, averaging the scores over multiple patients.

    Parameters:
    - df: The input DataFrame containing protein data.
    - list_of_proteins_to_predict: A list of proteins to predict.
    - proteins_list: The complete list of proteins in the DataFrame.
    - model_name: The name of the model to use for prediction ("XGBoost" or "Decision-Tree").
    - amount_of_patients: The number of patients to consider for averaging the scores.

    Returns:
    - cor_scores: A Counter object containing average correlation scores for each predicted protein.
    - r2_scores: A Counter object containing average R-squared scores for each predicted protein.
    - predictions: A Counter object containing averaged predictions for each protein.
    """
    df = df.copy()
    cor_scores, r2_scores, predictions = Counter(list()), Counter(list()), Counter(list())
    random_patients = random.sample(range(1, amount_of_patients + 1), min(5, amount_of_patients))
    for patient_number in random_patients:
        df_train = df.loc[
            df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
        df_test = df.loc[df['SampleID'] == patient_number]  # takes only patient patient_number for test

        for protein in list_of_proteins_to_predict:
            # predict one protein , we will put it inside Y_train:
            y_train, y_test = df_train[protein], df_test[protein]
            # we will put all the rest proteins inside X_train:
            pl_copy = proteins_list.copy()
            pl_copy.remove(protein)
            X_train, X_test = df_train[pl_copy], df_test[pl_copy]

            if (model_name == 'XGBoost'):
                cor_score, r2_score, prediction = model_XGBoostRegressor(X_train, y_train, X_test, y_test)
            else:
                cor_score, r2_score, prediction = model_DecisionTreeRegressor(X_train, y_train, X_test,
                                                                              y_test)

            cor_scores[protein] = cor_scores.get(protein, []) + [float(abs(cor_score[0, 1]))]
            if (r2_score < 0):
                r2_score = 0
            r2_scores[protein] = r2_scores.get(protein, []) + [r2_score]
            predictions[protein] = predictions.get(protein, []) + [prediction]

    for protein, value_list in cor_scores.items():
        cor_scores[protein] = sum(cor_scores[protein]) / len(cor_scores[protein])
        r2_scores[protein] = sum(r2_scores[protein]) / len(r2_scores[protein])
    return cor_scores, r2_scores, predictions


def model_DecisionTreeRegressor(X_train, y_train, X_test, y_test):
    """
    Trains a Decision Tree regressor model using the provided training data and makes predictions on the test data.

    Parameters:
    - X_train: The input features for training.
    - y_train: The target variable for training.
    - X_test: The input features for testing.
    - y_test: The target variable for testing.

    Returns:
    - cor: The correlation scores between y_test and the predictions.
    - r2: The R-squared score between y_test and the predictions.
    - prediction: The predicted values.
    """
    regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
    prediction = regressor.predict(X_test)
    if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
        cor = np.zeros((2, 2))
    else:
        cor = calculate_correlation(y_test, prediction)
    r2 = calculate_r2_score(y_test, prediction)
    return cor, r2, prediction


def model_XGBoostRegressor(X_train, y_train, X_test, y_test):
    """
   Trains an XGBoost regressor model using the provided training data and makes predictions on the test data.

   Parameters:
   - X_train: The input features for training.
   - y_train: The target variable for training.
   - X_test: The input features for testing.
   - y_test: The target variable for testing.

   Returns:
   - cor: The correlation scores between y_test and the predictions.
   - r2: The R-squared score between y_test and the predictions.
   - prediction: The predicted values.
    """
    regressor = xgb.XGBRegressor(random_state=0).fit(X_train, y_train)
    prediction = regressor.predict(X_test)
    if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
        cor = np.zeros((2, 2))
    else:
        cor = calculate_correlation(y_test, prediction)
    r2 = calculate_r2_score(y_test, prediction)
    return cor, r2, prediction


def calculate_correlation(y_test, prediction):
    """
    Calculates the correlation scores between y_test and the predictions.

    Parameters:
    - y_test: The actual target values.
    - prediction: The predicted values.

    Returns:
    - The correlation matrix.
    """
    return np.corrcoef(y_test.to_numpy(), prediction)


def calculate_r2_score(y_test, prediction):
    """
    Calculates the R-squared score between y_test and the predictions.

    Parameters:
    - y_test: The actual target values.
    - prediction: The predicted values.

    Returns:
    - The R-squared score.
    """
    return r2_score(y_test.to_numpy(), prediction)


def plot_graph_r2(DTR_r2_scores, model_name):
    """
    Plot a bar graph of R-squared scores for different proteins.

    Parameters:
    - DTR_r2_scores (dict): A dictionary containing R-squared scores for different proteins.
    - model_name (str): The name of the model.

    Returns:
    - fig: The matplotlib figure object.
    """
    data = DTR_r2_scores
    proteins = list(data.keys())
    scores = list(data.values())

    fig = plt.figure(figsize=(10, 7))

    # creating the bar plot
    plt.bar(proteins, scores, color='blue', width=0.4)
    plt.xticks(rotation=90, ha='right')
    plt.ylim(-1, 1)

    plt.xlabel("Proteins")
    plt.ylabel("r2 score")
    plt.title(model_name + " Scores")
    return fig


def plot_graph_cor(DTR_cor_scores, model_name):
    """
    Plot a bar graph of correlation scores for different proteins.

    Parameters:
    - DTR_cor_scores (dict): A dictionary containing correlation scores for different proteins.
    - model_name (str): The name of the model.

    Returns:
    - fig: The matplotlib figure object.
    """
    data = DTR_cor_scores
    proteins = list(data.keys())
    scores = list(data.values())

    fig = plt.figure(figsize=(10, 7))

    # creating the bar plot
    plt.bar(proteins, scores, color='red', width=0.4)
    plt.xticks(rotation=90, ha='right')
    plt.ylim(0, 1)

    plt.xlabel("Proteins")
    plt.ylabel("Correlation score")
    plt.title(model_name + " Scores")
    return fig


def prediction_matrix_creation(DTR_prediction, df, patient_number, cellLabel_image):
    """
    Create a prediction matrix based on the Deep Tissue Regression (DTR) predictions.

    Parameters:
    - DTR_prediction (list): A list of DTR predictions.
    - df (DataFrame): The input DataFrame containing the data.
    - patient_number (int): The number representing the test patient.
    - cellLabel_image (ndarray): An ndarray representing the cell labels in the image.

    Returns:
    - protein_prediction (ndarray): A 2D ndarray representing the prediction matrix.
    """
    df = df.copy()
    protein_prediction = np.zeros((np.size(cellLabel_image,0), np.size(cellLabel_image,1)))

    patient_numer_df = df.loc[df['SampleID'] == patient_number]  # takes only the test patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage']]
    protein_cellLabel_df['prediction'] = list(DTR_prediction)

    for index, row in protein_cellLabel_df.iterrows():
        protein_prediction[cellLabel_image == int(row['cellLabelInImage'])] = float(row['prediction'])
    return protein_prediction


def real_protein_matrix_creation(df, patient, protein, cellLabel_image):
    """
    Create a real protein matrix based on the patient's protein values.

    Parameters:
    - df (DataFrame): The input DataFrame containing the data.
    - patient (int): The number representing the patient.
    - protein (str): The name of the protein.
    - cellLabel_image (ndarray): An ndarray representing the cell labels in the image.

    Returns:
    - real_protein_matrix (ndarray): A 2D ndarray representing the real protein matrix.
    """
    df = df.copy()
    patient_numer_df = df.loc[df['SampleID'] == patient]  # takes only the patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage', protein]]
    real_protein_matrix = np.zeros((np.size(cellLabel_image,0), np.size(cellLabel_image,1)))

    for index, row in protein_cellLabel_df.iterrows():
        real_protein_matrix[cellLabel_image == int(row['cellLabelInImage'])] = float(row[protein])
    return real_protein_matrix


def visual_prediction(df, proteins_to_predict, patient):
    """
    Perform visual prediction analysis for a patient.

    Parameters:
    - df (DataFrame): The input DataFrame containing the data.
    - proteins_to_predict (list): A list of proteins to predict.
    - patient (int): The number representing the patient.
    """
    df = df.copy()
    flag = True
    # get from user cellLabel image:
    while flag:
        try:
            patient_labeled_cell_data = fd.askopenfilename()  # choose celldata of the patient
            cellLabel_image = Image.open(patient_labeled_cell_data)
            cellLabel_image = np.array(cellLabel_image)  # matrix of labeled cell data
            flag = False
        except:
            print("incoreect path to celldata.tiff of the testing patient")

    DTR_scores, DTR_r2_scores, DTR_prediction = ranking_model(df, patient, proteins_to_predict)
    for protein, protein_prediction in DTR_prediction.items():  # DTR_prediction is a dictionary
        print(f'starting protein : {protein}')
        prediction_matrix = prediction_matrix_creation(protein_prediction, df, patient, cellLabel_image)
        # prediction matrix to image:
        save_img(prediction_matrix, f'protein_prediction_{patient}_{protein}')
    return


def main(viewer, df, model_name, proteins_list, amount_of_patients):
    DTR_cor_scores, DTR_r2_scores, DTR_prediction = ranking_model_avg(df, proteins_list, proteins_list, model_name,
                                                                      amount_of_patients)
    # r2 plot
    plt2 = plot_graph_r2(dict(DTR_r2_scores.most_common()), model_name)
    plt2.savefig("plot_r2_ranking.png", dpi=170)
    napari_image1 = imread('plot_r2_ranking.png')  # Reads an image from file
    viewer.add_image(napari_image1,
                     name='plot_r2_ranking')  # Adds the image to the viewer and give the image layer a name
    # cor plot:
    plt1 = plot_graph_cor(dict(DTR_cor_scores.most_common()), model_name)
    plt1.savefig("plot_cor_ranking.png", dpi=170)
    napari_image2 = imread('plot_cor_ranking.png')  # Reads an image from file
    viewer.add_image(napari_image2,
                     name='plot_cor_ranking')  # Adds the image to the viewer and give the image layer a name


