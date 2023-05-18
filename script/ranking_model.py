import os
from collections import Counter
import pandas as pd
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score
import numpy as np
import xgboost as xgb
from PIL import Image


def ranking_model(df, patient_number, list_of_proteins_to_predict, proteins_list, model_name):
    df = df.copy()
    DTR_cor_scores, DTR_r2_scores, DTR_predictions = Counter(), Counter(), Counter()
    print(f'testing patient number :{patient_number}\n')
    df['SampleID'] = df['SampleID'].astype(str)
    df_train = df.loc[df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
    df_test = df.loc[df['SampleID'] == str(patient_number)]  # takes only patient patient_number for test

    for protein in list_of_proteins_to_predict:
        # predict one protein , we will put it inside Y_train:
        y_train = df_train[protein]
        y_test = df_test[protein]
        print(f'predicting protein: {protein}')
        # we will put all the rest proteins inside X_train:
        pl_copy = proteins_list.copy()
        pl_copy.remove(protein)
        X_train = df_train[pl_copy]
        X_test = df_test[pl_copy]
        if (model_name == 'XGBoost'):
            DTR_cor_score, DTR_r2_score, DTR_prediction = model_XGBoostRegressor(X_train, y_train, X_test, y_test)
        else:
            DTR_cor_score, DTR_r2_score, DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test, y_test)
        print(f'DTR r2 score: {DTR_r2_score}')
        print(f'DTR cor score: {DTR_cor_score[0, 1]}\n')

        DTR_cor_scores[protein] = float(DTR_cor_score[0, 1])
        DTR_r2_scores[protein] = DTR_r2_score
        DTR_predictions[protein] = DTR_prediction

    return DTR_cor_scores, DTR_r2_scores, DTR_predictions


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


def prediction_matrix_creation(DTR_prediction, df, patient_number, cellLabel_image):
    print(f'inside prediction_matrix_creation: DTR_prediction:\n{DTR_prediction}')
    df = df.copy()
    protein_prediction = np.zeros((np.size(cellLabel_image, 0), np.size(cellLabel_image, 1)))    df['SampleID'] = df['SampleID'].astype(str)
    patient_numer_df = df.loc[df['SampleID'] == patient_number]  # takes only the test patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage']]
    print(f'patient_numer_df: {patient_numer_df}')
    protein_cellLabel_df['prediction'] = list(DTR_prediction)
    print(f'inside prediction_matrix_creation: protein_cellLabel_df:\n{protein_cellLabel_df}')


    for index, row in protein_cellLabel_df.iterrows():
        protein_prediction[cellLabel_image == int(row['cellLabelInImage'])] = float(row['prediction'])

    return protein_prediction


def real_protein_matrix_creation(df, patient, protein, cellLabel_image):
    df = df.copy()
    patient_numer_df = df.loc[df['SampleID'] == patient]  # takes only the patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage', protein]]
    real_protein_matrix = np.zeros((np.size(cellLabel_image, 0), np.size(cellLabel_image, 1)))
    for index, row in protein_cellLabel_df.iterrows():
        real_protein_matrix[cellLabel_image == int(row['cellLabelInImage'])] = float(row[protein])
    return real_protein_matrix


def save_img(matrix, file_name):
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    # Append a number to the filename if it already exists
    file_number = 0
    image_filename = f'{file_name}.tiff'
    while os.path.exists(image_filename):
        file_number += 1
        image_filename = f'{file_name}{file_number}.tiff'
    print(image_filename)
    # save image
    new_im.save(image_filename)
    return image_filename
