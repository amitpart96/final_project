import threading
import time
import tkinter as tk
from idlelib.tooltip import _tooltip
from tkinter import filedialog as fd, filedialog
import PIL.Image
import pandas as pd
from PIL.Image import Image
from dask.array.chunk import view
from imageio import imread

import napari
import napari_hello
import os
from magicgui import magicgui
from enum import Enum
from napari.utils.notifications import show_info
import napari_hello
from napari_hello import create_csv
from napari_hello import ranking_model
from napari_hello import find_anomaly
from napari_hello import new_experiment
from magicgui.widgets import Select

import numpy as np
from sklearn.preprocessing import MinMaxScaler

class Options_Patients(Enum):  # todo: fill all the 41 patients
    patient1 = '1'
    patient2 = '2'
    patient3 = '3'
    patient4 = '4'
    patient5 = '5'
    patient6 = '6'
    patient7 = '7'
    patient8 = '8'
    patient9 = '9'
    patient10 = '10'


class Options_Proteins(Enum):
    Beta_catenin = 'Beta catenin'
    catenin = 'catenin'
    CD3 = 'CD3'
    CD4 = 'CD4'
    CD8 = 'CD8'
    CD11b = 'CD11b'
    CD11c = 'CD11c'
    CD16 = 'CD16'
    CD20 = 'CD20'
    CD31 = 'CD31'
    CD45 = 'CD45'
    CD45RO = 'CD45RO'
    CD56 = 'CD56'
    CD63 = 'CD63'
    CD68 = 'CD68'
    CD138 = 'CD138'
    CD209 = 'CD209'
    dsDNA = 'dsDNA'
    EGFR = 'EGFR'
    FoxP3 = 'FoxP3'
    H3K9ac = 'H3K9ac'
    H3K27me3 = 'H3K27me3'
    HLA_class_1 = 'HLA_class_1'
    HLA_DR = 'HLA-DR'
    IDO = 'IDO'
    Keratin6 = 'Keratin6'
    Keratin17 = 'Keratin17'
    Ki67 = 'Ki67'
    Lag3 = 'Lag3'
    MPO = 'MPO'
    p53 = 'p53'
    PD1 = 'PD1'
    PD_L1 = 'PD-L1'
    phospho_S6 = 'phospho-S6'
    SMA = 'SMA'
    Vimentin = 'Vimentin'

class Options_Proteins_New_Experiment(Enum):
    Keratin6 = 'Keratin6'
    FoxP3 = 'FoxP3'
    dsDNA = 'dsDNA'


viewer = napari.Viewer()
patient_number = None
protein = None
df = None
segmentation = None
prediction_matrix = None
real_protein_matrix = None
std_real = None
file_name_std = None
layer_std = None
protein_prediction_options_new_exeriment = []
list_of_proteins_to_train = []

@magicgui(chooseProteins=dict(widget_type='Select', choices=Options_Proteins),call_button='Predict Proteins')
def proteins_predict(chooseProteins):
    proteins_list = [protein.name for protein in chooseProteins]
    print(proteins_list)

    if (len(proteins_list) == 0):
        show_info("please select proteins")
        return
    if df is None:
        show_info("upload csv first")
        return
    if patient_number is None:
        show_info("choose patient number first")
        return
    show_info("starting to predict proteins")
    list_of_proteins_to_predict = proteins_list #["CD45", "dsDNA", "Vimentin"]
    ranking_model.predict_k_proteins(viewer, df, patient_number, list_of_proteins_to_predict)
    show_info('done predict proteins')
    return
# proteins_predict.show()

@magicgui(call_button='Upload cellTable')
def upload_csv():
    root = tk.Tk()
    root.withdraw()
    try:
        global df
        filename = fd.askopenfilename(title="open cellData csv")
        print(filename)
        df = pd.read_csv(filename)
        show_info(f'cellTable uploaded successfully')
        ranking_model_button.setVisible(True)
        patient_selection_button.setVisible(True)
    except:
        show_info("add path to cellData.csv in the code")
    return


@magicgui(call_button='Upload Segmentation')
def upload_segmentation():
    root = tk.Tk()
    root.withdraw()
    global segmentation
    seg = fd.askopenfilename()  # choose segmentation image of the patient
    # segmentation = Image.open(seg)
    napari_image = imread(seg)  # Reads an image from file
    viewer.add_image(napari_image, name="segmentation")  # Adds the image to the viewer and give the image layer a name
    show_info(f'segmentation uploaded successfully')
    return


def finish_create_csv():
    show_info('created cellTable successfully')
    create_CSV_button.setVisible(True)
    ranking_model_button.setVisible(True)
    patient_selection_button.setVisible(True)


def exepc_in_create_csv():
    create_CSV_button.setVisible(True)


@magicgui(call_button='Create cellTable')
def create_CSV():
    # create_CSV_button.setVisible(False)
    show_info("Processing cellTable creation")
    thread = threading.Thread(target=create_csv.main)
    thread.start()
    return


@magicgui(call_button='Rank Proteins')
def rankingg_model():
    if df is None:
        show_info("upload csv first")
        return
    show_info("starting to rank proteins")
    ranking_model.main(viewer, df, patient_number)
    show_info('done rank proteins')
    return



@magicgui(call_button='Find Anomaly')
def protein_selection(select_protein: Options_Proteins):
    # Do something with image and list of selected options
    global protein
    protein = select_protein.value
    show_info(f'{protein} is chosen')
    global prediction_matrix
    global real_protein_matrix
    global std_real
    global file_name_std
    global layer_std
    if df is None:
        show_info("upload csv first")
        return
    if patient_number is None:
        show_info("choose patient number first")
        return
    if protein is None:
        show_info("choose protein first")
        return
    show_info("starting to find anomaly")
    prediction_matrix, real_protein_matrix, std_real, file_name_std, layer_std = find_anomaly.main(viewer, df,
                                                                                                   patient_number,
                                                                                                   protein)
    show_info('done find anomaly')
    return


@magicgui(call_button='Select Patient')
def patient_selection(patient_selection: Options_Patients):
    # Do something with image and list of selected options
    global patient_number
    patient_number = int(patient_selection.value)
    show_info(f'patient {patient_number} is chosen')
    show_info(f'please upload patient {patient_number} proteins channels')
    root = tk.Tk()
    root.withdraw()
    list_img = fd.askopenfilenames(title="select the proteins channels of the patient")
    colors = list(napari.utils.colormaps.AVAILABLE_COLORMAPS)
    color = 0
    for img in list_img:
        channel_image = imread(img)  # Reads an image from file
        img_name = os.path.basename(img)
        img_name = img_name.removesuffix('.tiff') + " Patient" + str(patient_number)
        viewer.add_image(channel_image, name=img_name, colormap=colors[color],
                         visible=False)  # Adds the image to the viewer and give the image layer a name
        color += 1
        if (color >= len(colors)):
            color = 0
    show_info(f'images uploaded successfully')
    protein_selection_button.setVisible(True)
    k_proteins_predict_button.setVisible(True)
    change_std_button.setVisible(True)
    return


@magicgui(call_button='Upload Images')
def upload_images():
    root = tk.Tk()
    root.withdraw()
    list_img = fd.askopenfilenames()  # choose celldata of the patient
    for img in list_img:
        napari_image = imread(img)  # Reads an image from file
        img_name = os.path.basename(img)
        if patient_number is not None:
            img_name = "Patient" + str(patient_number) + " " + img_name
        viewer.add_image(napari_image, name=img_name)  # Adds the image to the viewer and give the image layer a name
    show_info(f'images uploaded successfully')
    return


@magicgui(
    call_button="change std",
    slider_float={"widget_type": "FloatSlider", 'max': 10},
)
def widget_demo(slider_float=2):
    print(slider_float)
    find_anomaly.update_difference_matrix_std(viewer, prediction_matrix, real_protein_matrix, std_real, slider_float,
                                              file_name_std, layer_std)
    return

@magicgui(call_button='Upload cellTable new experiment')
def upload_csv_new_experiment():
    root = tk.Tk()
    root.withdraw()
    try:
        global df_new_experiment
        global protein_prediction_options_new_exeriment
        global list_of_proteins_to_train
        filename = fd.askopenfilename(title="open cellData csv - new experiment", filetypes = (("CSV Files","*.csv"),))
        print(filename)
        df_new_experiment = pd.read_csv(filename)
        global df_new_experiment_normalized
        global df_normalized

        # copy the data
        df_new_experiment_normalized = df_new_experiment.copy()
        # apply normalization techniques
        columns = ['Na'] # todo - change the list of the columns to normalize - check if we want 2 separate lists for each table
        df_normalized = df.copy()
        for column in columns:
            #normalize cellTable new experiment
            df_new_experiment_normalized[column] = MinMaxScaler().fit_transform(np.array(df_new_experiment_normalized[column]).reshape(-1, 1))
            # normalize cellTable old experiment
            df_normalized[column] = MinMaxScaler().fit_transform(np.array(df_normalized[column]).reshape(-1, 1))
        print("1")
        protein_prediction_options_new_exeriment_tmp = set(df_normalized.columns) - set(df_new_experiment.columns)
        protein_prediction_options_new_exeriment = list(protein_prediction_options_new_exeriment_tmp)
        list_of_proteins_to_train =  np.intersect1d(df_normalized.columns, df_new_experiment_normalized.columns)
        print("2")
        print(protein_prediction_options_new_exeriment)
        print(f'train on proteins: {list_of_proteins_to_train}')
        #my_widget.changed.connect(update_choices)
        proteins_predict_new_experiment_button.setVisible(True)
        show_info(f'cellTable new experiment uploaded successfully')
        #predict_new_experiment.setVisible(True)
    except:
        show_info("add path to cellData.csv in the code")
    return

@magicgui(call_button='Select Patient New Experiment')
def patient_selection_new_experiment(patient_selection_new_experiment: Options_Patients):
    # Do something with image and list of selected options
    global patient_number_new_experiment
    patient_number_new_experiment = int(patient_selection_new_experiment.value)
    show_info(f'patient {patient_number_new_experiment} is chosen - new experiment')
    #show_info(f'please upload patient {patient_number} proteins channels')
    #root = tk.Tk()
    #root.withdraw()
    #list_img = fd.askopenfilenames(title="select the proteins channels of the patient")
    #colors = list(napari.utils.colormaps.AVAILABLE_COLORMAPS)
    #color = 0
    #for img in list_img:
     #   channel_image = imread(img)  # Reads an image from file
     #   img_name = os.path.basename(img)
     #   img_name = img_name.removesuffix('.tiff') + " Patient" + str(patient_number)
     #   viewer.add_image(channel_image, name=img_name, colormap=colors[color],
      #                   visible=False)  # Adds the image to the viewer and give the image layer a name
      #  color += 1
      #  if (color >= len(colors)):
      #      color = 0
    # show_info(f'images uploaded successfully')
    # protein_selection_button.setVisible(True)
    # k_proteins_predict_button.setVisible(True)
    # change_std_button.setVisible(True)
    return

@magicgui(choose_Proteins_New_Experiment=dict(widget_type='Select', choices=Options_Proteins_New_Experiment),call_button='Predict Proteins New Experiment')
def proteins_predict_new_experiment(choose_Proteins_New_Experiment):
    proteins_list = [protein.name for protein in choose_Proteins_New_Experiment]
    print(proteins_list)

    print(protein_prediction_options_new_exeriment)

    if (len(proteins_list) == 0):
        show_info("please select proteins - new Experiment")
        return
    if df_new_experiment_normalized is None:
        show_info("upload csv first")
        return
    if patient_number_new_experiment is None:
        show_info("choose patient number first (New Experiment)")
        return
   # show_info('done find anomaly')

    show_info("starting to predict proteins New Experiment")
    list_of_proteins_to_predict = proteins_list  # ["CD45", "dsDNA", "Vimentin"]
    new_experiment.predict_k_proteins(viewer, df_normalized,df_new_experiment_normalized,patient_number_new_experiment, list_of_proteins_to_predict, list_of_proteins_to_train)
    show_info('done predict proteins')
    return

#@magicgui(choices=["Item 1", "Item 2", "Item 3"])
#def my_widget(selected_items: ListWidget):
   # pass

#def update_choices(widget):
   # new_choices = protein_prediction_options_new_exeriment
  #  print(new_choices)
   # widget.update_choices(new_choices)



# widget_demo.show()
upload_segmentation_button = viewer.window.add_dock_widget(upload_segmentation, area='right')
create_CSV_button = viewer.window.add_dock_widget(create_CSV, area='right')
upload_csv_button = viewer.window.add_dock_widget(upload_csv, area='right')
ranking_model_button = viewer.window.add_dock_widget(rankingg_model, area='right')
patient_selection_button = viewer.window.add_dock_widget(patient_selection, area='right')
protein_selection_button = viewer.window.add_dock_widget(protein_selection, area='right')
change_std_button = viewer.window.add_dock_widget(widget_demo, area='right')
k_proteins_predict_button = viewer.window.add_dock_widget(proteins_predict, area='right')
upload_images_button = viewer.window.add_dock_widget(upload_images, area='right')
upload_csv_new_experiment_button = viewer.window.add_dock_widget(upload_csv_new_experiment, area='right')
patient_selection_new_experiment_button = viewer.window.add_dock_widget(patient_selection_new_experiment, area='right')
proteins_predict_new_experiment_button = viewer.window.add_dock_widget(proteins_predict_new_experiment, area='right')
#viewer.window.add_dock_widget(my_widget, area='right')

patient_selection_button.setVisible(False)
upload_images_button.setVisible(False)
ranking_model_button.setVisible(False)
k_proteins_predict_button.setVisible(False)
protein_selection_button.setVisible(False)
change_std_button.setVisible(False)
#upload_csv_new_experiment_button.setVisible(False)
#patient_selection_new_experiment_button.setVisible(False)
proteins_predict_new_experiment_button.setVisible(False)

def message():
    show_info('Welcome to Napari Plugin')


def main():
    napari


if __name__ == "__main__":
    main()
