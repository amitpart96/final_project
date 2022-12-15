import time
from collections import Counter
from random import randint
import numpy as np
import pandas as pd
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.metrics import r2_score


def ranking_model():
    GBR_scores, DTR_cor_scores, DTR_r2_scores = Counter(), Counter(), Counter()
    GBR_time, DTR_time = [], []
    df = pd.read_csv("napari_hello/cellData.csv")

    value = random_int = randint(1, 44)  # random chooses patient for test
    print(f'testing patient number {value}:\n')

    df_train = df.loc[df['SampleID'] != value]  # takes all patients for train, without patient value for test
    df_test = df.loc[df['SampleID'] == value]  # takes only patient value for test

    proteins_list = ["dsDNA", "Vimentin", "SMA", "FoxP3", "Lag3", "CD4", "CD16", "CD56", "PD1", "CD31", "PD-L1", "EGFR",
                     "Ki67", "CD209", "CD11c", "CD138", "CD68", "CD8", "CD3", "Keratin17", "IDO", "CD63", "CD45RO", "CD20",
                     "p53", "Beta catenin", "HLA-DR", "CD11b", "CD45", "H3K9ac", "Pan-Keratin", "H3K27me3",
                     "phospho-S6", "MPO", "Keratin6", "HLA_Class_1"]

    for protein in proteins_list:
        # predict one protein , we will put it inside Y_train:
        y_train, y_test = df_train[protein], df_test[protein]
        print(f'predicting protein: {protein}')
        # we will put all the rest proteins inside X_train:
        pl_copy = proteins_list.copy()
        pl_copy.remove(protein)
        X_train, X_test = df_train[pl_copy], df_test[pl_copy]

        # DecisionTreeRegressor:
        startDTR = time.time()
        DTR_cor_score, DTR_r2_score, DTR_prediction = model_DecisionTreeRegressor(X_train, y_train, X_test, y_test)
        finishDTR = time.time()
        DTR_time.append(finishDTR - startDTR)
        print(f'DTR r2 score: {DTR_r2_score}')
        print(f'DTR cor score: {DTR_cor_score[0, 1]}\n')
        #print("DTR prediction: " + str(DTR_prediction))


        '''
        # GradientBoostingRegressor:
        startGBR = time.time()
        GBR_score, GBR_prediction = model_GradientBoostingRegressor(X_train, y_train, X_test, y_test)
        finishGBR = time.time()
        GBR_time.append(finishGBR - startGBR)
        print("GBR score: " + str(GBR_score[0, 1]))
        # print("GBR prediction: " + str(GBR_prediction))
        GBR_scores[protein] = GBR_score[0, 1]
        '''
        DTR_cor_scores[protein] = float(DTR_cor_score[0, 1])
        DTR_r2_scores[protein] = DTR_r2_score

    return DTR_cor_scores, DTR_r2_scores, sum(DTR_time)


def model_GradientBoostingRegressor(X_train, y_train, X_test, y_test):
    est = GradientBoostingRegressor(n_estimators=100, learning_rate=0.1, max_depth=1, random_state=0,
                                    loss='squared_error').fit(X_train, y_train)
    prediction = est.predict(X_test)
    if np.std(y_test.to_numpy()) == 0 or np.std(prediction) == 0:
        print(
            "The correlation could not be computed because the standard deviation of one of the series is equal to zero")
        cor = np.zeros((2, 2))
    else:
        cor = calculate_correlation(y_test, prediction)
    r2 = calculate_r2_score(y_test, prediction)
    return cor, r2, prediction


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


if __name__ == "__main__":
    start = time.time()
    DTR_scores, DTR_r2_scores, DTR_time = ranking_model()
    # ranked_proteins_GBR = sorted(GBR_scores, key=GBR_scores.get, reverse=True)
    ranked_proteins_DTR_by_cor = sorted(DTR_scores, key=DTR_scores.get, reverse=True)
    ranked_proteins_DTR_by_r2 = sorted(DTR_r2_scores, key=DTR_r2_scores.get, reverse=True)

    # print(f'ranked_proteins_GBR:\n{ranked_proteins_GBR}')
    # print(f'total time of GBR: {GBR_time}\n')

    print(f'ranked_proteins_DTR_by_cor:\n{ranked_proteins_DTR_by_cor}')
    print(f'ranked_proteins_DTR_by_r2:\n{ranked_proteins_DTR_by_r2}')
    # print(f'total time of DTR: {DTR_time}\n')

    end = time.time()
    total_time = end - start

    print(f'total time: {total_time}')
    # V להוסיף מדד r2 score (ולהדרג לפיו) מ-sklearn
    # לקחת 3 פציאנטים ולעשו תעליהם ממוצע
    # להכניס לפלאגאין
    # V ניתאי ישלח רשימה מעודכנת של חלבונים
    #


# Press the green button in the gutter to run the script.
def main():
    start = time.time()
    DTR_scores, DTR_r2_scores, DTR_time = ranking_model()
    # ranked_proteins_GBR = sorted(GBR_scores, key=GBR_scores.get, reverse=True)
    ranked_proteins_DTR_by_cor = sorted(DTR_scores, key=DTR_scores.get, reverse=True)
    ranked_proteins_DTR_by_r2 = sorted(DTR_r2_scores, key=DTR_r2_scores.get, reverse=True)

    # print(f'ranked_proteins_GBR:\n{ranked_proteins_GBR}')
    # print(f'total time of GBR: {GBR_time}\n')

    print(f'ranked_proteins_DTR_by_cor:\n{ranked_proteins_DTR_by_cor}')
    print(f'ranked_proteins_DTR_by_r2:\n{ranked_proteins_DTR_by_r2}')
    #print(f'total time of DTR: {DTR_time}\n')

    end = time.time()
    total_time = end - start

    print(f'total time: {total_time}')
    # V להוסיף מדד r2 score (ולהדרג לפיו) מ-sklearn
    # לקחת 3 פציאנטים ולעשו תעליהם ממוצע
    # להכניס לפלאגאין
    # V ניתאי ישלח רשימה מעודכנת של חלבונים
    #
