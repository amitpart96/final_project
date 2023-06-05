from sklearn.tree import DecisionTreeRegressor
import pandas as pd
from tkinter import filedialog as fd
import numpy as np
from PIL import Image
from skimage.io import imread
import xgboost as xgb

def predict_k_proteins(viewer, df,df_new_experiment,patient_number_new_experiment, list_of_proteins_to_predict, list_of_proteins_to_train,patient_cellLabel_image, model_name):
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
    # we will put all the rest proteins inside X_train:
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
        #for pycharm run test, uncomment the next 2 rows:
        # cellLabel_image = Image.open(img)
        # print(np.asarray(cellLabel_image))
        #comment the next row when checking in napari
        viewer.add_image(protein_prediction_image, name=img_name)  # Adds the image to the viewer and give the image layer a name
    return DTR_prediction

# normilizing the series
def min_max_scaling(series):
    if series.max() - series.min() !=0:
        return (series - series.min()) / (series.max() - series.min())
    return pd.Series(0, index=range(series.size))

#Predict for DecisionTree model
def model_DecisionTreeRegressor(X_train, y_train, X_test):
    regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
    prediction = regressor.predict(X_test)
    return prediction

#Predict for XGBoost model
def model_XGBoostRegressor(X_train, y_train, X_test):
    regressor = xgb.XGBRegressor(random_state=0).fit(X_train, y_train)  # Instantiate XGBoost regressor
    prediction = regressor.predict(X_test)
    return prediction


def find_anomaly(viewer,df, df_new_experiment,patient_number_new_experiment,protein_to_predict,list_of_proteins_to_train,patient_cellLabel_image):  # todo: change name to find anomaly
    df_normalized = df.copy()
    df_new_experiment_normalized = df_new_experiment.copy()
    print(f'starting patient number: {patient}')
    scores, r2_scores, prediction = predict(df_normalized,df_new_experiment_normalized, patient_number_new_experiment, protein_to_predict,list_of_proteins_to_train, model_name)

    for protein, protein_prediction in prediction.items():  # DTR_prediction is a dictionary
        print(f'starting protein : {protein}')
        prediction_matrix = ranking_model.prediction_matrix_creation(protein_prediction, df, patient, patient_cellLabel_image)
        # prediction matrix to image:
        pred_img = ranking_model.save_img(prediction_matrix, f'protein_prediction_{patient}_{protein}')
        real_protein_matrix = ranking_model.real_protein_matrix_creation(df, patient, protein, patient_cellLabel_image)
        # real matrix to image:
        real_img = ranking_model.save_img(real_protein_matrix, f'real_protein_{patient}_{protein}')
        # std of the real matrix
        std_real = real_protein_matrix.std()
        print(std_real)
        # difference by std
        difference_matrix_std = create_difference_matrix_std(prediction_matrix, real_protein_matrix, std_real)
        difference_matrix = abs(np.subtract(real_protein_matrix, prediction_matrix))
        # difference_matrix to image:
        diff_img = ranking_model.save_img(difference_matrix, f'difference_matrix_{patient}_{protein}')
        diff_img_std = ranking_model.save_img(difference_matrix_std, f'difference_matrix_std_{patient}_{protein}')
        file_name_std = f'difference_matrix_std_{patient}_{protein}'
    print(f'finished patient number: {patient}')
    return real_img, pred_img, diff_img, diff_img_std, prediction_matrix, real_protein_matrix, std_real, file_name_std  # ??

def predict(df_normalized,df_new_experiment_normalized, patient_number_new_experiment, protein_to_predict,list_of_proteins_to_train, model_name):
    df_normalized = df_normalized.copy()
    df_new_experiment_normalized = df_new_experiment_normalized.copy()
    DTR_cor_scores, DTR_r2_scores, DTR_predictions = Counter(), Counter(), Counter()

    print(f'testing patient number :{patient_number}\n')
    df_train = df.loc[df['SampleID'] != patient_number]  # takes all patients for train, without patient patient_number for test
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




# def model_DecisionTreeRegressor(X_train, y_train, X_test):
#     regressor = DecisionTreeRegressor(random_state=0).fit(X_train, y_train)
#     # cross_val_score(regressor, X_train, y_train, cv=10)
#     prediction = regressor.predict(X_test)
#     #if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
#       #  print(
#        #     "The correlation could not be computed because the standard deviation of one of the series is equal to zero")
#        # cor = np.zeros((2, 2))
#     #else:
#       #  cor = calculate_correlation(y_test, prediction)
#    # r2 = calculate_r2_score(y_test, prediction)
#     return prediction

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

def save_img(matrix, file_name):
    matrix = (matrix * 255).round().astype(np.uint8)
    new_im = Image.fromarray(matrix)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    print(image_filename)
    # save image using extension
    new_im.save(image_filename)
    return image_filename