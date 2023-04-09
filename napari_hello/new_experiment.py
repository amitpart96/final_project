from sklearn.tree import DecisionTreeRegressor


def predict_k_proteins(viewer, df_normalized,df_new_experiment_normalized,patient_number_new_experiment, list_of_proteins_to_predict, list_of_proteins_to_train):
    df_normalized = df_normalized.copy()
    df_new_experiment_normalized = df_new_experiment_normalized.copy()
    print(f'testing patient number :{patient_number_new_experiment}\n')
    df_train = df_normalized.loc[df_normalized['SampleID'] != patient_number_new_experiment]  # takes all patients for train, from the df_normalized
    df_test = df_new_experiment_normalized.loc[df_new_experiment_normalized['SampleID'] == patient_number_new_experiment]  # takes only patient patient_number_new_experiment for test from df_new_experiment_normalized

    proteins_list = ["CD45", "dsDNA", "Vimentin", "SMA", "FoxP3", "Lag3", "CD4", "CD16", "CD56", "PD1", "CD31",
                     "PD-L1",
                     "EGFR",
                     "Ki67", "CD209", "CD11c", "CD138", "CD68", "CD8", "CD3", "Keratin17", "IDO", "CD63", "CD45RO",
                     "CD20",
                     "p53", "Beta catenin", "HLA-DR", "CD11b", "H3K9ac", "Pan-Keratin", "H3K27me3",
                     "phospho-S6", "MPO", "Keratin6", "HLA_Class_1"]

    list_of_proteins_to_train = ['Au', 'C','CD16'  ,'CD20', 'CD209', 'CD3', 'CD31', 'CD4'
 ,'CD45' ,'CD68', 'CD8' ,'IDO' ,'Ki67' ,'Lag3' ,'MPO' ,'Na' ,'SMA', 'SampleID' ,'Ta'
 ,'Vimentin'] # 'CD11b' ,'CD11c' ,'CD163'
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
    return DTR_prediction

'''
    flag = True
    # get from user cellLabel image:
    while flag:
        try:
            patient_labeled_cell_data = fd.askopenfilename(title=f'choose cellData image of the patient {patient_number}')  # choose celldata of the patient
            cellLabel_image = Image.open(patient_labeled_cell_data)
            cellLabel_image = np.array(cellLabel_image)  # matrix of labeled cell data
            flag = False
        except:
            print("incoreect path to celldata.tiff of the testing patient")
    prediction_df = pd.DataFrame(DTR_prediction, columns=list_of_proteins_to_predict)

    for protein_name, values in prediction_df.iteritems():
        protein_prediction_matrix = prediction_matrix_creation(prediction_df[protein_name], df, patient_number, cellLabel_image)
        print(f'protein_prediction_matrix:\n {protein_prediction_matrix}')

        img_name= f'protein_prediction_{patient_number}_{protein_name}'
        img = save_img(protein_prediction_matrix, img_name)
        protein_prediction_image = imread(img)
        #for pycharm run test, uncomment the next 2 rows:
        # cellLabel_image = Image.open(img)
        # print(np.asarray(cellLabel_image))
        #comment the next row when checking in napari
        viewer.add_image(protein_prediction_image, name=img_name)  # Adds the image to the viewer and give the image layer a name
    return DTR_cor_score, DTR_r2_score, DTR_prediction


'''



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