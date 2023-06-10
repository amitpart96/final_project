import threading
import time
import tkinter as tk
from idlelib.tooltip import _tooltip
from tkinter import filedialog as fd, filedialog
import pandas as pd
from PIL import Image
from dask.array.chunk import view
from imageio import imread
import napari
import napari_hello
import os
from magicgui import magicgui
from magicgui import magic_factory
from enum import Enum
from napari.utils.notifications import show_info
import napari_hello
from napari_hello import create_csv
from napari_hello import ranking_model
from napari_hello import find_anomaly
from napari_hello import new_experiment
from napari_hello import predict_k_proteins
from napari_hello import new_experiment
from napari_hello import segmentation
from magicgui.widgets import Select
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import re
from napari.utils.notifications import show_warning
from pathlib import Path
from napari.types import LayerDataTuple


class Options_Models(Enum):
    DecisionTree = 'DecisionTree'
    XGBoost = 'XGBoost'


viewer = napari.Viewer()
patient_number = None
protein = None
df = None
prediction_matrix = None
real_protein_matrix = None
std_real = None
file_name_std = None
layer_std = None
protein_prediction_options_new_exeriment = []
list_of_proteins_to_train = []
model_name = None
cellLabelImages = None
patient_cellLabel_image = None
root_directory_path = None

DEFAULT_CHOICES = []
DEFAULT_CHOICES_find_anomaly = []
DEFAULT_CHOICES_PATIENT = []
func_flag = False

def _update_choices_on_file_change(widget):
    """Return available choices of proteins from csv headers"""
    filename = widget.filename.value
    choices = None
    choices_patient = None
    global df_new_experiment
    widget.dropdown.choices = choices = ()
    widget.dropdown_patient.choices = ()
    if filename.is_file():
        if filename.suffix == ".csv":
            df_new_experiment = pd.read_csv(filename)
            choices = sorted(set(df.columns) - set(df_new_experiment.columns))
            choices_patient = sorted(df_new_experiment["SampleID"].drop_duplicates())
        else:
            show_warning(f"File {filename} is not a .csv file.")
    if choices is not None:
        widget.dropdown.choices = choices
    if choices_patient is not None:
        #print(type(widget.dropdown_patient.choices))
        #print(widget.dropdown_patient.choices)
        #widget.dropdown_patient.choices = choices_patient
        #print(tuple(sorted(widget.dropdown_patient.choices)))
        #widget.dropdown_patient.choices = tuple(sorted(widget.dropdown_patient.choices))
        widget.dropdown_patient.choices = choices_patient
        print(widget.dropdown_patient.choices)

def proteins_predict1():
    @magicgui(
        filename={
            "label": "CSV file with proteins to predict",
            "mode": "r",
            "filter": "*.csv",
        },
        foldername={
            "label": "folder with cellLabeled Images of the patients",
            "mode": "d",
        },
        dropdown_patient=dict(
            widget_type="ComboBox", choices=DEFAULT_CHOICES_PATIENT, label="Patient to predict",
        ),
        dropdown=dict(
            widget_type="Select", choices=DEFAULT_CHOICES, label="Proteins to predict",
        ),
        comboBox_model=dict(
            widget_type="ComboBox", choices=['DesicionTree', 'XGBoost'], label="model Selection",
        ),
        call_button="Predict Proteins New Experiment",
        
    )

    def widget(viewer: napari.Viewer, dropdown_patient, dropdown,comboBox_model,filename=Path.home(), foldername=Path.home()):
        # Perform the prediction here
        if df is None:
            show_info("upload csv first")
            return
        if df_new_experiment is None:
            show_info("upload csv new experiment")
            return
        if (len(dropdown) == 0):
            show_info("please select proteins - new Experiment")
            return

        if dropdown_patient is None:
            show_info("choose patient number first (New Experiment)")
            return



        proteins_list_to_predict = dropdown
        patient_number = dropdown_patient
        print(proteins_list_to_predict)
        print(patient_number)
        patient_number_new_experiment = patient_number
        #proteins_list = [protein.name for protein in choose_Proteins_New_Experiment]

        # list of proteins to train - the intersection of the proteins of old experiment and new experiment
        list_of_proteins_to_train = np.setdiff1d(np.intersect1d(df.columns, df_new_experiment.columns),np.array(['SampleID', 'cellLabelInImage', 'cellSize']))

   # show_info('done find anomaly')
        print(foldername)
        # loop through all the files in the selected directory

        image_paths = []
        for filename in os.listdir(foldername):
            # check if the file is an image file
            if filename.endswith(".tiff") | filename.endswith(".tif"):
                image_path = os.path.join(foldername, filename)
                image_paths.append(image_path)
        print("before find patient_cellLabel_image")
        print(image_paths)
        print(patient_number_new_experiment)
        patient_cellLabel_image_new_experiment = find_cell_labeled_image(image_paths, patient_number_new_experiment)

        if patient_cellLabel_image_new_experiment is None:
            show_info("choose folder with cellLabeled images of the patient (New Experiment)")
            return

        show_info("starting to predict proteins New Experiment")
        list_of_proteins_to_predict = proteins_list_to_predict  # ["CD45", "dsDNA", "Vimentin"]
        print(patient_number_new_experiment)
        print(proteins_list_to_predict)
        print(list_of_proteins_to_train)
        new_experiment.predict_k_proteins(viewer, df,df_new_experiment,patient_number_new_experiment, proteins_list_to_predict, list_of_proteins_to_train, patient_cellLabel_image_new_experiment, comboBox_model)
        _update_choices_on_file_change(widget)
        show_info('done predict proteins')



    @widget.filename.changed.connect
    def update_choices_on_file_change(event=None):
        _update_choices_on_file_change(widget)


    return widget


# @magicgui(call_button='Upload Segmentation')
# def upload_segmentation():
#     root = tk.Tk()
#     root.withdraw()
#     global segmentation
#     seg = fd.askopenfilename()  # choose segmentation image of the patient
#     # segmentation = Image.open(seg)
#     napari_image = imread(seg)  # Reads an image from file
#     viewer.add_image(napari_image, name="segmentation")  # Adds the image to the viewer and give the image layer a name
#     show_info(f'segmentation uploaded successfully')
#     return


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
    thread = threading.Thread(target=create_csv.main(root_directory_path))
    thread.start()
    return


@magicgui(call_button='Choose Model')
def choose_model(model_selection: Options_Models):
    global model_name
    model_name = model_selection.value
    show_info(f'{model_name} model is chosen')
    return


def get_proteins_list(df):
    column_names = df.columns.tolist()
    # Columns to remove
    columns_to_remove = ['SampleID', 'cellLabelInImage', 'cellSize']  # List of column names to remove

    # Remove columns from the list
    list_of_proteins_from_df = [col for col in column_names if col not in columns_to_remove]

    return list_of_proteins_from_df


@magicgui(call_button='Rank Proteins')
def rankingg_model():
    if df is None:
        show_info("upload csv first")
        return
    if model_name is None:
        show_info("choose model first")
        return
    dir_contents = os.listdir(root_directory_path)
    amount_of_patients = len(dir_contents)
    show_info("starting to rank proteins")
    proteins_list = get_proteins_list(df)

    ranking_model.main(viewer, df, model_name, proteins_list, amount_of_patients)
    show_info('done rank proteins')
    return


@magicgui(call_button='Create Segmentation')
def create_seggmentation():
    show_info("starting to create segmentations")
    global root_directory_path
    print("root2 : ",root_directory_path)
    segmentation.main(root_directory_path)
    show_info('done to create segmentations')
    return


@magicgui(call_button='New Experiment')
def new_exp():
    create_segmentation_button.setVisible(True)
    create_CSV_button.setVisible(True)

    # upload_segmentation_button.setVisible(False)
    upload_csv_button.setVisible(False)
    patient_selection_button.setVisible(False)
    # upload_images_button.setVisible(False)
    choose_model_button.setVisible(False)
    ranking_model_button.setVisible(False)
    k_proteins_predict_button.setVisible(False)
    protein_selection_button.setVisible(False)
    other_exp_m_button.setVisible(False)

    new_exp_button.setVisible(False)
    old_exp_button.setVisible(True)
    return


@magicgui(call_button='Old Experiment')
def old_exp():
    # upload_segmentation_button.setVisible(True)
    upload_csv_button.setVisible(True)

    create_segmentation_button.setVisible(False)
    create_CSV_button.setVisible(False)
    new_exp_button.setVisible(True)
    old_exp_button.setVisible(False)
    return


def upload_exp():
    @magicgui(
        foldername={
            "label": "root experiment directory",
            "mode": "d",
        },
        call_button="Upload New Experiment",
    )
    # todo:
    def widget(viewer: napari.Viewer, foldername=Path.home()):
        if Path.home() != foldername:
            global root_directory_path
            # open the directory dialog box and allow the user to select a directory
            root_directory_path = foldername
            show_info(f'root experiment directory chosen successfully')
            upload_exp_button.setVisible(False)
            new_exp_button.setVisible(True)
            old_exp_button.setVisible(True)

    return widget


def find_cell_labeled_image(paths, patient_number):
    for path in paths:
        file_name = os.path.basename(path)

        match_1 = re.search(r"p(\d+)_labeledcellData\.tiff", file_name)
        match_2 = re.search(r"p(\d+)_labeledcellData\.tif", file_name)
        if match_1:
            img_number_1 = int(match_1.group(1))
            if img_number_1 == patient_number:
                print(patient_number)
                img = Image.open(path)
                array_img = np.asarray(img)
                return array_img
        elif match_2:
            img_number_2 = int(match_2.group(1))
            if img_number_2 == patient_number:
                print(patient_number)
                img = Image.open(path)
                array_img = np.asarray(img)
                return array_img
    return None


exp_proteins_choices = None
exp_patients_choices = None
DEFAULT_CHOICES_PATIENTS = []


def patient_selection():
    @magicgui(
        dropdown=dict(
            widget_type="ComboBox", choices=DEFAULT_CHOICES_PATIENTS, label="Patient:",
        ),
        call_button='Select Patient',
    )
    def widget(viewer: napari.Viewer, dropdown):
        if df is None:
            print("return patient widget")
            return widget

        global patient_number
        patient_number = int(dropdown)
        show_info(f'patient {patient_number} is chosen')

        # uploading to the napari viewer the proteins images of the patient
        directory_path = os.path.join(root_directory_path, str(patient_number))
        # Get a list of all files and directories in the directory
        files = os.listdir(directory_path)
        list_img = [(os.path.join(directory_path, file)) for file in files]
        colors = list(napari.utils.colormaps.AVAILABLE_COLORMAPS)
        color = 0
        for img in list_img:
            channel_image = imread(img)  # Reads an image from file
            img_name = os.path.basename(img)
            img_name = img_name + " Patient" + str(patient_number)
            viewer.add_image(channel_image, name=img_name, colormap=colors[color],
                             visible=False)
            color += 1
            if color >= len(colors):
                color = 0
        show_info(f'patient images uploaded successfully')
        global patient_cellLabel_image
        patient_cellLabel_image = find_cell_labeled_image(cellLabelImages, patient_number)
        protein_selection_button.setVisible(True)
        k_proteins_predict_button.setVisible(True)

    return widget


DEFAULT_CHOICES_find_anomalies = []


def find_anomlayy():
    @magicgui(
        dropdown=dict(
            widget_type="ComboBox", choices=DEFAULT_CHOICES_find_anomalies, label="Protein to find anomalies",
        ),
        slider_float={"widget_type": "FloatSlider", 'max': 10},
        call_button="Find anomalies",
    )
    def widget(viewer: napari.Viewer, dropdown, slider_float=2):
        if df is None:
            show_info("upload csv first")
            return widget
        if cellLabelImages is None:
            show_info("upload cellLabelImages first")
            return widget
        if patient_number is None:
            show_info("choose patient number first")
            return widget
        if model_name is None:
            show_info("choose model first")
            return widget
        if patient_cellLabel_image is None:
            show_info("couldn't find the cellLabelImagePath, make sure the image name is correct")
            return widget
        global protein
        protein = dropdown
        show_info(f'{protein} is chosen')
        global prediction_matrix
        global real_protein_matrix
        global std_real
        global file_name_std
        global layer_std
        show_info("starting to find anomaly")
        proteins_list = get_proteins_list(df)
        prediction_matrix, real_protein_matrix, std_real, file_name_std, layer_std = \
            find_anomaly.main(viewer, df, patient_number, protein, model_name, proteins_list, patient_cellLabel_image, slider_float)

    @widget.slider_float.changed.connect #new
    def update_choices_on_slider_float_change(event=None): #new
        _update_choices_on_slider_float_change(widget) #new

    return widget

def _update_choices_on_slider_float_change(widget):
    """Return available choices of proteins from csv headers"""
    slider_float = widget.slider_float.value
    try:
        find_anomaly.update_difference_matrix_std(viewer, prediction_matrix, real_protein_matrix, std_real, slider_float, file_name_std, layer_std)
    except:
        return






def upload_CellTable_and_cellLabelImage(protein_widget, patients_widget, find_anomaly_widget):
    @magicgui(
        filename={
            "label": "CSV file with proteins to predict",
            "mode": "r",
            "filter": "*.csv",
        },
        foldername={
            "label": "folder with cellLabeled Images of the patients",
            "mode": "d",
        },
        call_button="Update files Uploading",
    )
    # todo:
    def widget(viewer: napari.Viewer, filename=Path.home(), foldername=Path.home()):
        if (Path.home() != filename):
            global df
            df = pd.read_csv(filename)
            show_info(f'cellTable uploaded successfully')
        if (Path.home() != foldername):
            global cellLabelImages
            cellLabelImages = []
            for filename in os.listdir(foldername):
                image_path = os.path.join(foldername, filename)
                cellLabelImages.append(image_path)
            show_info(f'cellLabelImages uploaded successfully')
        if (filename != Path.home() and foldername != Path.home()):
            choose_model_button.setVisible(True)
            ranking_model_button.setVisible(True)
            patient_selection_button.setVisible(True)
            other_exp_m_button.setVisible(True)
            _update_proteins_choices_on_file_changes(widget, protein_widget)
            _update_proteins_choices_on_file_changes(widget, find_anomaly_widget)
            _update_patient_choices_on_file_changes(widget, patients_widget)

    return widget


DEFAULT_CHOICES_OF_PROTEINS = []


def _update_proteins_choices_on_file_changes(filename_widget, dropdown_widget):
    """Return available choices of proteins from csv headers"""
    filename = filename_widget.filename.value
    choices = None
    if filename.is_file():
        if filename.suffix == ".csv":
            df = pd.read_csv(filename)
            choices = get_proteins_list(df)
        else:
            show_warning(f"File {filename} is not a .csv file.")
    if choices is not None:
        dropdown_widget.dropdown.choices = choices
        global exp_proteins_choices
        exp_proteins_choices = choices


def _update_patient_choices_on_file_changes(filename_widget, patients_widget):
    """Return available choices of patients from csv """
    filename = filename_widget.filename.value
    choices = None
    if filename.is_file():
        if filename.suffix == ".csv":
            df = pd.read_csv(filename)
            choices = df["SampleID"].drop_duplicates()
        else:
            show_warning(f"File {filename} is not a .csv file.")
    if choices is not None:
        patients_widget.dropdown.choices = choices
        global exp_patients_choices
        exp_patients_choices = choices
    return


def _update_proteins_choices_on_file_changes2(df, widget):
    choices = get_proteins_list(df)
    if choices is not None:
        widget.dropdown.choices = choices


def _update_patient_choices_on_file_changes2(df, widget):
    choices = df["SampleID"].drop_duplicates()
    if choices is not None:
        widget.dropdown.choices = choices


def proteins_predict():
    @magicgui(
        dropdown=dict(
            widget_type="Select", choices=DEFAULT_CHOICES_OF_PROTEINS, label="Proteins to predict"
        ),
        call_button="Predict Proteins",
    )
    def widget(viewer: napari.Viewer, dropdown):
        # Perform the prediction here
        proteins_list_to_predict = dropdown

        if (len(proteins_list_to_predict) == 0):
            show_info("please select proteins")
            return widget
        if df is None:
            show_info("upload csv first")
            return widget
        if patient_cellLabel_image is None:
            show_info("couldn't find the cellLabelImagePath, make sure the image name is correct")
            return widget
        if patient_number is None:
            show_info("choose patient number first")
            return widget
        if model_name is None:
            show_info("choose model first")
            return widget
        show_info("starting to predict proteins")
        list_of_proteins_to_predict = proteins_list_to_predict
        proteins_list = get_proteins_list(df)
        predict_k_proteins.predict_k_proteins(viewer, df, patient_number, list_of_proteins_to_predict, proteins_list,
                                              model_name, patient_cellLabel_image)
        show_info('done predict proteins')

    return widget


protein_widget = proteins_predict()
find_anomaly_widget = find_anomlayy()
patients_widget = patient_selection()

# upload, new, old experiment buttons:
upload_exp_button = viewer.window.add_dock_widget(upload_exp(), area='right', name="Upload New Experiment")
new_exp_button = viewer.window.add_dock_widget(new_exp, area='right', name="New Experiment")
old_exp_button = viewer.window.add_dock_widget(old_exp, area='right', name="Old Experiment")

# create buttons:
create_segmentation_button = viewer.window.add_dock_widget(create_seggmentation, area='right',
                                                           name="Create Segmentation")
create_CSV_button = viewer.window.add_dock_widget(create_CSV, area='right', name="Create new CellTable")

# uplode buttons:
# upload_segmentation_button = viewer.window.add_dock_widget(upload_segmentation, area='right')
upload_csv_button = viewer.window.add_dock_widget(
    upload_CellTable_and_cellLabelImage(protein_widget, patients_widget, find_anomaly_widget),
    name="Upload cellTable and cellLabelImages")

choose_model_button = viewer.window.add_dock_widget(choose_model, area='right', name="Choose Model")
ranking_model_button = viewer.window.add_dock_widget(rankingg_model, area='right', name="Rank Proteins")
patient_selection_button = viewer.window.add_dock_widget(patients_widget, area='right', name="Patient Selection")
patient_selection_button.setToolTip('Click on upload cellTable and cellLabelImages to refresh options')
protein_selection_button = viewer.window.add_dock_widget(find_anomaly_widget, area='right', name="Find Anomalies")
protein_selection_button.setToolTip('Click on upload cellTable and cellLabelImages to refresh options')

k_proteins_predict_button = viewer.window.add_dock_widget(protein_widget, name="Predict Proteins", area='right')

other_exp_m_button = viewer.window.add_dock_widget(proteins_predict1(), name="Predict proteins new experiment",
                                                   area='right')

new_exp_button.setVisible(False)
old_exp_button.setVisible(False)
create_segmentation_button.setVisible(False)
create_CSV_button.setVisible(False)
# upload_segmentation_button.setVisible(False)
upload_csv_button.setVisible(False)
patient_selection_button.setVisible(False)
choose_model_button.setVisible(False)
ranking_model_button.setVisible(False)
k_proteins_predict_button.setVisible(False)
protein_selection_button.setVisible(False)
other_exp_m_button.setVisible(False)


def message():
    show_info('Welcome to Napari Plugin')


def main():
    napari


if __name__ == "__main__":
    main()
