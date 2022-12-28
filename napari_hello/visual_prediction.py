from random import randint
from turtle import fd
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
from PIL.Image import Image
from napari_hello import ranking_model
from collections import Counter
from random import randint
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image


def visual_prediction(df, proteins_to_predict, patient):
    df = df.copy()
    print(f'starting patient number: {patient}')
    flag=True
    while flag:
        try:
            print("entered try")
            patient_labeled_cell_data = fd.askopenfilename()  # choose celldata of the patient
            print("1")
            cellLabel_image = Image.open(patient_labeled_cell_data)
            cellLabel_image = np.array(cellLabel_image)
            print("2")
            #cellLabel_image = np.array(cellLabel_image)  # matrix of labeled cell data
            print("3")
            flag = False
            print("done try")
        except:
            print("incoreect path to celldata.tiff of the testing patient")

    print("before ranking model function call")
    DTR_scores, DTR_r2_scores, DTR_prediction = ranking_model.ranking_model(df, patient, proteins_to_predict)
    print("after ranking model function call")

    for protein, protein_prediction in DTR_prediction.items():  # DTR_prediction is a dictionary
        print(f'starting protein : {protein}')
        prediction_matrix = ranking_model.prediction_matrix_creation(protein_prediction, df, patient, cellLabel_image)
        # prediction matrix to image:

        ranking_model.save_img(prediction_matrix, f'protein_prediction_{patient}_{protein}')
    print(f'finished patient number: {patient}')


def main():
    root = tk.Tk()
    root.withdraw()

    try:
        filename = fd.askopenfilename()
        print(filename)
        df = pd.read_csv(filename)

    except:
        print("add path to cellData.csv in the code")

    patient_number = random_int = randint(1, 42)  # random chooses patient for test
    #list_of_proteins_to_predict = top5 = ['CD45', 'CD45RO', 'CD3', 'CD4', 'H3K27me3']
    list_of_proteins_to_predict = top5 = ['CD45']
    visual_prediction(df, list_of_proteins_to_predict,patient_number)

if __name__ == "__main__":
    main()
