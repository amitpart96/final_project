from collections import Counter
from random import randint
import pandas as pd
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
import matplotlib.pyplot as plt
from skimage.io import imread

def save_img(matrix, file_name):
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    print(image_filename)
    # save image using extension
    new_im.save(image_filename)
    return image_filename

def ranking_model(df, patient_number, list_of_proteins_to_predict):
    df = df.copy()
    DTR_cor_scores, DTR_r2_scores, DTR_predictions = Counter(), Counter(), Counter()

    print(f'testing patient number :{patient_number}\n')
    df_train = df.loc[
        df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
    df_test = df.loc[df['SampleID'] == patient_number]  # takes only patient patient_number for test

    proteins_list = ["CD45", "dsDNA", "Vimentin", "SMA", "FoxP3", "Lag3", "CD4", "CD16", "CD56", "PD1", "CD31", "PD-L1",
                     "EGFR",
                     "Ki67", "CD209", "CD11c", "CD138", "CD68", "CD8", "CD3", "Keratin17", "IDO", "CD63", "CD45RO",
                     "CD20",
                     "p53", "Beta catenin", "HLA-DR", "CD11b", "H3K9ac", "Pan-Keratin", "H3K27me3",
                     "phospho-S6", "MPO", "Keratin6", "HLA_Class_1"]

    for protein in list_of_proteins_to_predict:  # need to change to top 5 list
        # predict one protein , we will put it inside Y_train:
        y_train, y_test = df_train[protein], df_test[protein]
        print(f'predicting protein: {protein}')
        # we will put all the rest proteins inside X_train:
        pl_copy = proteins_list.copy()
        pl_copy.remove(protein)
        X_train, X_test = df_train[pl_copy], df_test[pl_copy]

        # DecisionTreeRegressor:
        DTR_cor_score, DTR_r2_score, DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test, y_test)
        print(f'DTR r2 score: {DTR_r2_score}')
        print(f'DTR cor score: {DTR_cor_score[0, 1]}\n')
        # print("DTR prediction: " + str(DTR_prediction))

        DTR_cor_scores[protein] = float(DTR_cor_score[0, 1])
        DTR_r2_scores[protein] = DTR_r2_score
        DTR_predictions[protein] = DTR_prediction
    return DTR_cor_scores, DTR_r2_scores, DTR_predictions


def model_DecisionTreeRegressor(X_train, y_train, X_test, y_test):
    regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
    # cross_val_score(regressor, X_train, y_train, cv=10)
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
    return np.corrcoef(y_test.to_numpy(), prediction)


def calculate_r2_score(y_test, prediction):
    return r2_score(y_test.to_numpy(), prediction)


def calculate_r2_score(y_test, prediction):
    return r2_score(y_test.to_numpy(), prediction)


def plot_graph_r2(DTR_r2_scores):
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
    plt.title("Decision Tree Regressor Scores")
    #plt.show()
    return plt


def plot_graph_cor(DTR_cor_scores):
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
    plt.title("Decision Tree Regressor Scores")
    #plt.show()
    return fig


def prediction_matrix_creation(DTR_prediction, df, patient_number, cellLabel_image):
    df = df.copy()
    protein_prediction = np.zeros((2048, 2048))

    patient_numer_df = df.loc[df['SampleID'] == patient_number]  # takes only the test patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage']]
    protein_cellLabel_df["prediction"] = DTR_prediction

    for index, row in protein_cellLabel_df.iterrows():
        protein_prediction[cellLabel_image == int(row['cellLabelInImage'])] = float(row['prediction'])

    return protein_prediction


def real_protein_matrix_creation(df, patient, protein, cellLabel_image):
    df = df.copy()
    patient_numer_df = df.loc[df['SampleID'] == patient]  # takes only the patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage', protein]]
    real_protein_matrix = np.zeros((2048, 2048))

    for index, row in protein_cellLabel_df.iterrows():
        real_protein_matrix[cellLabel_image == int(row['cellLabelInImage'])] = float(row[protein])
    return real_protein_matrix


def five_patients_prediction(df, top_5_proteins):
    df = df.copy()
    for patient in range(1, 6):  # need to chaznge to random 5 patients ???
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

        DTR_scores, DTR_r2_scores, DTR_prediction = ranking_model(df, patient, top_5_proteins)
        for protein, protein_prediction in DTR_prediction.items():  # DTR_prediction is a dictionary
            print(f'starting protein : {protein}')
            prediction_matrix = prediction_matrix_creation(protein_prediction, df, patient, cellLabel_image)
            # prediction matrix to image:
            save_img(prediction_matrix, f'protein_prediction_{patient}_{protein}')
            # real protein matrix creation:
            real_protein_matrix = real_protein_matrix_creation(df, patient, protein, cellLabel_image)
            # real matrix to image:
            save_img(real_protein_matrix, f'real_protein_{patient}_{protein}')
            difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
            # difference_matrix to image:
            save_img(difference_matrix, f'difference_matrix_{patient}_{protein}')
            # save all 3 images for each protein and each patient ???
            print(f'finished 3 images for protein: {protein}')
        print(f'finished patient number: {patient}')
    return


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


def find_the_best_pro(df, protein_list):
    df = df.copy()
    array_r2_scores = []
    for patient in range(1, 16):
        print(f'starting patient number: {patient}')
        DTR_scores, DTR_r2_scores, DTR_prediction = ranking_model(df, patient, protein_list)
        ranked_proteins_DTR_by_r2 = sorted(DTR_r2_scores, key=DTR_r2_scores.get, reverse=True)[:7]
        array_r2_scores.append(ranked_proteins_DTR_by_r2)
    print(f'array_r2_scores:{array_r2_scores}')


def main(viewer,df,patient_number):
    proteins_list = ["CD45", "dsDNA", "Vimentin"]
  # proteins_list = ["CD45", "dsDNA", "Vimentin", "SMA", "FoxP3", "Lag3", "CD4", "CD16", "CD56", "PD1", "CD31", "PD-L1",
  #                    "EGFR",
  #                    "Ki67", "CD209", "CD11c", "CD138", "CD68", "CD8", "CD3", "Keratin17", "IDO", "CD63", "CD45RO",
  #                    "CD20",
  #                    "p53", "Beta catenin", "HLA-DR", "CD11b", "H3K9ac", "Pan-Keratin", "H3K27me3",
  #                    "phospho-S6", "MPO", "Keratin6", "HLA_Class_1"]


    DTR_cor_scores, DTR_r2_scores, DTR_prediction = ranking_model(df, patient_number, proteins_list)

    #r2 plot
    plt2 = plot_graph_r2(dict(DTR_r2_scores.most_common()))
    plt2.savefig("plot_r2_ranking.png", dpi=150)
    napari_image1 = imread('plot_r2_ranking.png')  # Reads an image from file
    viewer.add_image(napari_image1, name='plot_r2_ranking')  # Adds the image to the viewer and give the image layer a name
    #cor plot:
    plt1 = plot_graph_cor(dict(DTR_cor_scores.most_common()))
    plt1.savefig("plot_cor_ranking.png", dpi=150)
    napari_image2 = imread('plot_cor_ranking.png')  # Reads an image from file
    viewer.add_image(napari_image2, name='plot_cor_ranking')  # Adds the image to the viewer and give the image layer a name

    # ranked_proteins_DTR_by_cor = sorted(DTR_scores, key=DTR_scores.get, reverse=True)
    # ranked_proteins_DTR_by_r2 = sorted(DTR_r2_scores, key=DTR_r2_scores.get, reverse=True)
    # print(f'ranked_proteins_DTR_by_cor:\n{ranked_proteins_DTR_by_cor}')
    # print(f'ranked_proteins_DTR_by_r2:\n{ranked_proteins_DTR_by_r2}')


if __name__ == "__main__":
    main()
