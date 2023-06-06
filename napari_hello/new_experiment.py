from sklearn.tree import DecisionTreeRegressor
import pandas as pd
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
from skimage.io import imread
import xgboost as xgb

def predict_k_proteins(viewer, df,df_new_experiment,patient_number_new_experiment, list_of_proteins_to_predict, list_of_proteins_to_train,patient_cellLabel_image, model_name):
    """
    Predict protein expression levels for a set of proteins using a trained model.

    Parameters:
    - viewer: The Napari viewer object.
    - df (pd.DataFrame): The training data DataFrame containing protein expression levels for training.
    - df_new_experiment (pd.DataFrame): The experimental data DataFrame containing protein expression levels for prediction.
    - patient_number_new_experiment: The patient number of the new experiment for prediction.
    - list_of_proteins_to_predict (list): A list of protein names to predict their expression levels.
    - list_of_proteins_to_train (list): A list of protein names used for training the prediction model.
    - patient_cellLabel_image: The path to the cellLabel image of the testing patient.
    - model_name (str): The name of the prediction model to use ('XGBoost' or 'DecisionTree').

    Returns:
    np.ndarray: An array of predicted protein expression levels.
    """
    df = df.copy()
    df_new_experiment = df_new_experiment.copy()
    #normalization of the proteins in the experiment to values 0-1
    for col in list(list_of_proteins_to_train)+list_of_proteins_to_predict:
        df[col] = min_max_scaling(df[col])
    for col in list_of_proteins_to_train:
        df_new_experiment[col] = min_max_scaling(df_new_experiment[col])

    df_normalized = df
    df_new_experiment_normalized = df_new_experiment

    df_train = df_normalized # takes all patients for train, from the df_normalized
    df_test = df_new_experiment_normalized.loc[df_new_experiment_normalized['SampleID'] == patient_number_new_experiment]  # takes only patient patient_number_new_experiment for test from df_new_experiment_normalized

    # predict one protein , we will put it inside Y_train:
    if (len(list_of_proteins_to_predict) == 1):
        y_train = df_train[list_of_proteins_to_predict[0]]
    else: # y_train for multiple proteins prediction
        y_train = df_train[list_of_proteins_to_predict]
    pl_copy = list_of_proteins_to_train.copy()
    X_train, X_test = df_train[pl_copy], df_test[pl_copy]


    if (model_name == 'XGBoost'):
        DTR_prediction = model_XGBoostRegressor(X_train, y_train, X_test) #predicion in case of XGBoost model
    else:
        DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test) #predicion in case of DecisionTree model

    print(f'DTR prediction: {DTR_prediction}')
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
        viewer.add_image(protein_prediction_image, name=img_name)  # Adds the image to the viewer and give the image layer a name
    return DTR_prediction

def min_max_scaling(series):
    """
    Perform min-max scaling on a pandas Series.

    Parameters:
    - series (pd.Series): The input Series to be scaled.

    Returns:
    pd.Series: The scaled Series.
    """
    if series.max() - series.min() !=0:
        return (series - series.min()) / (series.max() - series.min())
    return pd.Series(0, index=range(series.size))

def model_DecisionTreeRegressor(X_train, y_train, X_test):
    """
    Perform regression using the Decision Tree Regressor model.

    Parameters:
    - X_train (pd.DataFrame): The training data features.
    - y_train (pd.Series): The training data target values.
    - X_test (pd.DataFrame): The test data features.

    Returns:
    np.ndarray: The predicted target values for the test data.
    """
    regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
    prediction = regressor.predict(X_test)
    return prediction

def model_XGBoostRegressor(X_train, y_train, X_test):
    """
    Perform regression using the XGBoost Regressor model.

    Parameters:
    - X_train (pd.DataFrame): The training data features.
    - y_train (pd.Series): The training data target values.
    - X_test (pd.DataFrame): The test data features.

    Returns:
    np.ndarray: The predicted target values for the test data.
    """
    regressor = xgb.XGBRegressor(random_state=0).fit(X_train, y_train)  # Instantiate XGBoost regressor
    prediction = regressor.predict(X_test)
    return prediction

def prediction_matrix_creation(DTR_prediction, df, patient_number, cellLabel_image):
    """
    Create a prediction matrix based on the predicted values and cell labels.

    Parameters:
    - DTR_prediction (np.ndarray): The predicted values for the proteins.
    - df (pd.DataFrame): The dataframe containing information about the cells.
    - patient_number (str or int): The patient number for which predictions are being made.
    - cellLabel_image (np.ndarray): The image matrix containing cell labels.

    Returns:
    np.ndarray: The prediction matrix where each cell is assigned a predicted protein value.
    """
    df = df.copy()
    protein_prediction = np.zeros((np.size(cellLabel_image,0), np.size(cellLabel_image,1)))

    patient_numer_df = df.loc[df['SampleID'] == patient_number]  # takes only the test patient
    protein_cellLabel_df = patient_numer_df[['cellLabelInImage']]
    protein_cellLabel_df['prediction'] = list(DTR_prediction)

    for index, row in protein_cellLabel_df.iterrows():
        protein_prediction[cellLabel_image == int(row['cellLabelInImage'])] = float(row['prediction'])

    return protein_prediction

def save_img(matrix, file_name):
    """
    Save a matrix as an image file.

    Parameters:
    - matrix (np.ndarray): The matrix to be saved as an image.
    - file_name (str): The name of the image file to be saved.

    Returns:
    str: The filename of the saved image.
    """
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    # save image using extension
    new_im.save(image_filename)
    return image_filename