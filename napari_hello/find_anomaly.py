from turtle import fd
from PIL.Image import Image
from imageio import imread
from napari_hello import ranking_model
import pandas as pd
import tkinter as tk
from tkinter import filedialog as fd
import numpy as np
from PIL import Image


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
            cellLabel_image = np.array(cellLabel_image)
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
        difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
        # difference_matrix to image:
        diff_img=ranking_model.save_img(difference_matrix, f'difference_matrix_{patient}_{protein}')
    print(f'finished patient number: {patient}')
    return real_img,pred_img,diff_img

def main(viewer):
    root = tk.Tk()
    root.withdraw()

    try:
        filename = fd.askopenfilename()
        print(filename)
        df = pd.read_csv(filename)

    except:
        print("add path to cellData.csv in the code")

    #patient_number = random_int = randint(1, 42)  # random chooses patient for test
    patient_number = 1
    #list_of_proteins_to_predict = top5 = ['CD45', 'CD45RO', 'CD3', 'CD4', 'H3K27me3']
    list_of_proteins_to_predict = top5 = ['CD45']
    real_img,pred_img,diff_img=find_anomaly(df, list_of_proteins_to_predict, patient_number)
    napari_image = imread(real_img)  # Reads an image from file
    viewer.add_image(napari_image, name=real_img)  # Adds the image to the viewer and give the image layer a name
    napari_image = imread(pred_img)  # Reads an image from file
    viewer.add_image(napari_image, name=pred_img)  # Adds the image to the viewer and give the image layer a name
    napari_image = imread(diff_img)  # Reads an image from file
    viewer.add_image(napari_image, name=diff_img)  # Adds the image to the viewer and give the image layer a name

if __name__ == "__main__":
    main()
