from sklearn.tree import DecisionTreeRegressor
import pandas as pd
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
from skimage.io import imread

def predict_k_proteins(viewer, df_normalized,df_new_experiment_normalized,patient_number_new_experiment, list_of_proteins_to_predict, list_of_proteins_to_train,patient_cellLabel_image):
    df_normalized = df_normalized.copy()
    df_new_experiment_normalized = df_new_experiment_normalized.copy()
    print(f'testing patient number :{patient_number_new_experiment}\n')
    df_train = df_normalized.loc[df_normalized['SampleID'] != patient_number_new_experiment]  # takes all patients for train, from the df_normalized
    df_test = df_new_experiment_normalized.loc[df_new_experiment_normalized['SampleID'] == patient_number_new_experiment]  # takes only patient patient_number_new_experiment for test from df_new_experiment_normalized

    '''proteins_list = ["CD45", "dsDNA", "Vimentin", "SMA", "FoxP3", "Lag3", "CD4", "CD16", "CD56", "PD1", "CD31",
                     "PD-L1",
                     "EGFR",
                     "Ki67", "CD209", "CD11c", "CD138", "CD68", "CD8", "CD3", "Keratin17", "IDO", "CD63", "CD45RO",
                     "CD20",
                     "p53", "Beta catenin", "HLA-DR", "CD11b", "H3K9ac", "Pan-Keratin", "H3K27me3",
                     "phospho-S6", "MPO", "Keratin6", "HLA_Class_1"]'''

    # predict one protein , we will put it inside Y_train:
    if (len(list_of_proteins_to_predict) == 1):
        y_train = df_train[list_of_proteins_to_predict[0]]
    else:
        y_train = df_train[list_of_proteins_to_predict]
    # we will put all the rest proteins inside X_train:
    pl_copy = list_of_proteins_to_train.copy()
    #for protein in list_of_proteins_to_predict:
       # pl_copy.remove(protein)
    X_train, X_test = df_train[pl_copy], df_test[pl_copy]

    # DecisionTreeRegressor:
    DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test)
    print(f'DTR prediction: {DTR_prediction}')
    #todo: לשאול את ניתאי איך להציג את הקורלציה של יותר מחלבון אחד
    #print(f'DTR cor score: {DTR_cor_score}\n')
    #with open(..., 'w', newline='') as myfile:
      #  wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
       # wr.writerow(mylist)
    #return DTR_prediction
    try:
        cellLabel_image = patient_cellLabel_image
    except:
        print("incoreect path to celldata.tiff of the testing patient")

    prediction_df = pd.DataFrame(DTR_prediction, columns=list_of_proteins_to_predict)

    for protein_name, values in prediction_df.iteritems():
        protein_prediction_matrix = prediction_matrix_creation(prediction_df[protein_name], df_new_experiment_normalized, patient_number_new_experiment, cellLabel_image)
        print(f'protein_prediction_matrix:\n {protein_prediction_matrix}')

        img_name= f'protein_prediction_{patient_number_new_experiment}_{protein_name}'
        img = save_img(protein_prediction_matrix, img_name)
        protein_prediction_image = imread(img)
        #for pycharm run test, uncomment the next 2 rows:
        # cellLabel_image = Image.open(img)
        # print(np.asarray(cellLabel_image))
        #comment the next row when checking in napari
        viewer.add_image(protein_prediction_image, name=img_name)  # Adds the image to the viewer and give the image layer a name
    return DTR_prediction






def model_DecisionTreeRegressor(X_train, y_train, X_test):
    regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
    # cross_val_score(regressor, X_train, y_train, cv=10)
    prediction = regressor.predict(X_test)
    #if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
      #  print(
       #     "The correlation could not be computed because the standard deviation of one of the series is equal to zero")
       # cor = np.zeros((2, 2))
    #else:
      #  cor = calculate_correlation(y_test, prediction)
   # r2 = calculate_r2_score(y_test, prediction)
    return prediction

def prediction_matrix_creation(DTR_prediction, df, patient_number, cellLabel_image):
    print(f'inside prediction_matrix_creation: DTR_prediction:\n{DTR_prediction}')
    df = df.copy()
    protein_prediction = np.zeros((1024, 1024))

    patient_numer_df = df.loc[df['SampleID'] == patient_number]  # takes only the test patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage']]
    protein_cellLabel_df['prediction'] = list(DTR_prediction)
    print(f'inside prediction_matrix_creation: protein_cellLabel_df:\n{protein_cellLabel_df}')

    for index, row in protein_cellLabel_df.iterrows():
        protein_prediction[cellLabel_image == int(row['cellLabelInImage'])] = float(row['prediction'])

    return protein_prediction

def save_img(matrix, file_name):
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    print(image_filename)
    # save image using extension
    new_im.save(image_filename)
    return image_filename