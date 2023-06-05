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
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    print(image_filename)
    # save image using extension
    new_im.save(image_filename)
    return image_filename


def ranking_model(df, patient_number, list_of_proteins_to_predict, proteins_list, model_name):
    df = df.copy()
    DTR_cor_scores, DTR_r2_scores, DTR_predictions = Counter(), Counter(), Counter()

    print(f'testing patient number :{patient_number}\n')
    df_train = df.loc[
        df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
    df_test = df.loc[df['SampleID'] == patient_number]  # takes only patient patient_number for test

    for protein in list_of_proteins_to_predict:
        # predict one protein , we will put it inside Y_train:
        y_train, y_test = df_train[protein], df_test[protein]
        print(f'predicting protein: {protein}')
        # we will put all the rest proteins inside X_train:
        pl_copy = proteins_list.copy()
        pl_copy.remove(protein)
        X_train, X_test = df_train[pl_copy], df_test[pl_copy]

        if (model_name == 'XGBoost'):
            DTR_cor_score, DTR_r2_score, DTR_prediction = model_XGBoostRegressor(X_train, y_train, X_test, y_test)
        else:
            DTR_cor_score, DTR_r2_score, DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test, y_test)
        print(f'DTR r2 score: {DTR_r2_score}')
        print(f'DTR cor score: {DTR_cor_score[0, 1]}\n')
        # print("DTR prediction: " + str(DTR_prediction))

        DTR_cor_scores[protein] = float(DTR_cor_score[0, 1])
        DTR_r2_scores[protein] = DTR_r2_score
        DTR_predictions[protein] = DTR_prediction

    return DTR_cor_scores, DTR_r2_scores, DTR_predictions


def ranking_model_avg(df, list_of_proteins_to_predict, proteins_list, model_name, amount_of_patients):
    print(model_name)
    df = df.copy()
    cor_scores, r2_scores, predictions = Counter(list()), Counter(list()), Counter(list())
    random_patients = random.sample(range(1, amount_of_patients + 1), min(5, amount_of_patients))
    for patient_number in random_patients:
        print(f'testing patient number :{patient_number}\n')
        df_train = df.loc[
            df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
        df_test = df.loc[df['SampleID'] == patient_number]  # takes only patient patient_number for test

        for protein in list_of_proteins_to_predict:
            # predict one protein , we will put it inside Y_train:
            y_train, y_test = df_train[protein], df_test[protein]
            print(f'predicting protein: {protein}')
            # we will put all the rest proteins inside X_train:
            pl_copy = proteins_list.copy()
            pl_copy.remove(protein)
            X_train, X_test = df_train[pl_copy], df_test[pl_copy]

            if (model_name == 'XGBoost'):
                cor_score, r2_score, prediction = model_XGBoostRegressor(X_train, y_train, X_test, y_test)
            else:
                cor_score, r2_score, prediction = model_DecisionTreeRegressor(X_train, y_train, X_test,
                                                                              y_test)
            print(f'r2 score: {r2_score}')
            print(f'cor score: {cor_score[0, 1]}\n')
            # print("DTR prediction: " + str(DTR_prediction))

            cor_scores[protein] = cor_scores.get(protein, []) + [float(abs(cor_score[0, 1]))]
            if (r2_score < 0):
                r2_score = 0
            r2_scores[protein] = r2_scores.get(protein, []) + [r2_score]
            predictions[protein] = predictions.get(protein, []) + [prediction]
    print(f'cor_scores dict before avg: {cor_scores}')
    print(f'r2_scores dict before avg: {r2_scores}')

    for protein, value_list in cor_scores.items():
        cor_scores[protein] = sum(cor_scores[protein]) / len(cor_scores[protein])
        r2_scores[protein] = sum(r2_scores[protein]) / len(r2_scores[protein])
    print(f'cor_scores dict after avg: {cor_scores}')
    print(f'r2_scores dict after avg: {r2_scores}')
    return cor_scores, r2_scores, predictions


def model_DecisionTreeRegressor(X_train, y_train, X_test, y_test):
    regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
    prediction = regressor.predict(X_test)
    if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
        print(
            "The correlation could not be computed because the standard deviation of one of the series is equal to zero")
        cor = np.zeros((2, 2))
    else:
        cor = calculate_correlation(y_test, prediction)
    r2 = calculate_r2_score(y_test, prediction)
    return cor, r2, prediction


def model_XGBoostRegressor(X_train, y_train, X_test, y_test):
    regressor = xgb.XGBRegressor(random_state=0).fit(X_train, y_train)  # Instantiate XGBoost regressor
    prediction = regressor.predict(X_test)
    if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
        print(
            "The correlation could not be computed because the standard deviation of one of the series is equal to zero")
        cor = np.zeros((2, 2))
    else:
        cor = calculate_correlation(y_test, prediction)
    r2 = calculate_r2_score(y_test, prediction)
    return cor, r2, prediction


def calculate_correlation(y_test, prediction):
    print(f'y_test.to_numpy():{y_test.to_numpy()}')
    print(f'prediction: {prediction}')
    return np.corrcoef(y_test.to_numpy(), prediction)


def calculate_r2_score(y_test, prediction):
    return r2_score(y_test.to_numpy(), prediction)


def plot_graph_r2(DTR_r2_scores, model_name):
    # creating the dataset
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
    # plt.show()
    return fig


def plot_graph_cor(DTR_cor_scores, model_name):
    # creating the dataset
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
    # plt.show()
    return fig


def prediction_matrix_creation(DTR_prediction, df, patient_number, cellLabel_image):
    print(f'inside prediction_matrix_creation: DTR_prediction:\n{DTR_prediction}')
    df = df.copy()
    protein_prediction = np.zeros((np.size(cellLabel_image,0), np.size(cellLabel_image,1)))

    patient_numer_df = df.loc[df['SampleID'] == patient_number]  # takes only the test patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage']]
    protein_cellLabel_df['prediction'] = list(DTR_prediction)
    print(f'inside prediction_matrix_creation: protein_cellLabel_df:\n{protein_cellLabel_df}')

    for index, row in protein_cellLabel_df.iterrows():
        protein_prediction[cellLabel_image == int(row['cellLabelInImage'])] = float(row['prediction'])

    return protein_prediction


def real_protein_matrix_creation(df, patient, protein, cellLabel_image):
    df = df.copy()
    patient_numer_df = df.loc[df['SampleID'] == patient]  # takes only the patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage', protein]]
    real_protein_matrix = np.zeros((np.size(cellLabel_image,0), np.size(cellLabel_image,1)))

    for index, row in protein_cellLabel_df.iterrows():
        real_protein_matrix[cellLabel_image == int(row['cellLabelInImage'])] = float(row[protein])
    return real_protein_matrix


def visual_prediction(df, proteins_to_predict, patient):
    df = df.copy()
    print(f'starting patient number: {patient}')
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
    print(f'finished patient number: {patient}')
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


if __name__ == "__main__":
