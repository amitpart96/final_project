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
# from napari_hello import segmentation
from magicgui.widgets import Select
import numpy as np
from sklearn.preprocessing import MinMaxScaler
import re
from napari.utils.notifications import show_warning
from pathlib import Path

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


class Options_Models(Enum):
    DecisionTree = 'DecisionTree'
    XGBoost = 'XGBoost'

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
model_name = None
cellLabelImages = None
patient_cellLabel_image = None
root_directory_path = None


@magicgui(chooseProteins=dict(widget_type='Select', choices=Options_Proteins), call_button='Predict Proteins')
def proteins_predict(chooseProteins):
    proteins_list_to_predict = [protein.name for protein in chooseProteins]
    print(proteins_list_to_predict)

    if (len(proteins_list_to_predict) == 0):
        show_info("please select proteins")
        return
    if df is None:
        show_info("upload csv first")
        return
    if patient_cellLabel_image is None:
        show_info("couldn't find the cellLabelImagePath, make sure the image name is correct")
        return
    if patient_number is None:
        show_info("choose patient number first")
        return
    if model_name is None:
        show_info("choose model first")
        return
    show_info("starting to predict proteins")
    list_of_proteins_to_predict = proteins_list_to_predict  # ["CD45", "dsDNA", "Vimentin"]
    proteins_list = get_proteins_list(df)
    predict_k_proteins.predict_k_proteins(viewer, df, patient_number, list_of_proteins_to_predict, proteins_list, model_name,patient_cellLabel_image)
    show_info('done predict proteins')
    return


# proteins_predict.show()

@magicgui(call_button='Upload cellTable and cellLabelImages')
def upload_csv():
    root = tk.Tk()
    root.withdraw()
    try:
        global df
        filename = fd.askopenfilename(title="open cellData csv")
        print(filename)
        df = pd.read_csv(filename)
        show_info(f'cellTable uploaded successfully')

    except:
        show_info("add path to cellData.csv in the code")

    root = tk.Tk()
    root.withdraw()

    global cellLabelImages
    try:
        # create an empty list to store the image data
        cellLabelImages = []
        # open the directory dialog box and allow the user to select a directory
        dir_path = filedialog.askdirectory(title="choose cellLabelImages directory")

        # loop through all the files in the selected directory
        for filename in os.listdir(dir_path):
            # check if the file is an image file
            if filename.endswith(".tiff"):
                image_path = os.path.join(dir_path, filename)
                cellLabelImages.append(image_path)
        print(f'cellLabelImages={cellLabelImages}')
        show_info(f'cellLabelImages uploaded successfully')
    except:
        show_info("add path to cellLabelImages directory")

    choose_model_button.setVisible(True)
    ranking_model_button.setVisible(True)
    patient_selection_button.setVisible(True)
    upload_csv_new_experiment_button.setVisible(True)
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
    if cellLabelImages is None:
        show_info("upload cellLabelImages first")
        return
    if patient_number is None:
        show_info("choose patient number first")
        return
    if protein is None:
        show_info("choose protein first")
        return
    if model_name is None:
        show_info("choose model first")
        return
    if patient_cellLabel_image is None:
        show_info("couldn't find the cellLabelImagePath, make sure the image name is correct")
        return
    show_info("starting to find anomaly")
    proteins_list = get_proteins_list(df)
    prediction_matrix, real_protein_matrix, std_real, file_name_std, layer_std = find_anomaly.main(viewer, df,
                                                                                                   patient_number,
                                                                                                   protein, model_name,
                                                                                                   proteins_list,
                                                                                                   patient_cellLabel_image)
    show_info('done find anomaly')
    return


@magicgui(call_button='Select Patient')
def patient_selection(patient_selection: Options_Patients):
    # Do something with image and list of selected options
    global patient_number
    patient_number = int(patient_selection.value)
    show_info(f'patient {patient_number} is chosen')

    #uploading to the napari viwer the proteins images of the patient
    directory_path = os.path.join(root_directory_path, str(patient_number))

    # Get a list of all files and directories in the directory
    files = os.listdir(directory_path)
    list_img = [(os.path.join(directory_path, file)) for file in files]

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
    show_info(f'patient images uploaded successfully')

    global patient_cellLabel_image
    patient_cellLabel_image = find_path_of_patient(cellLabelImages, patient_number)

    protein_selection_button.setVisible(True)
    k_proteins_predict_button.setVisible(True)
    change_std_button.setVisible(True)
    return


def find_path_of_patient(paths, patient_number):
    for path in paths:
        file_name = os.path.basename(path)

        match = re.search(r"p(\d+)_labeledcellData\.tiff", file_name)

        if match:
            img_number = int(match.group(1))
            if img_number == patient_number:
                img = Image.open(path)
                array_img = np.asarray(img)
                return array_img
    return None




@magicgui(
    call_button="change std",
    slider_float={"widget_type": "FloatSlider", 'max': 10},
)
def widget_demo(slider_float=2):
    print(slider_float)
    find_anomaly.update_difference_matrix_std(viewer, prediction_matrix, real_protein_matrix, std_real, slider_float,
                                              file_name_std, layer_std)
    return


#new


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

        show_info(f'cellTable new experiment uploaded successfully')

        patient_selection_new_experiment_button.setVisible(True)

    except:
        show_info("add path to cellData.csv in the code")
    return

@magicgui(call_button='Select Patient New Experiment')
def patient_selection_new_experiment(patient_selection_new_experiment: Options_Patients):
    # Do something with image and list of selected options
    global patient_number_new_experiment
    patient_number_new_experiment = int(patient_selection_new_experiment.value)
    show_info(f'patient {patient_number_new_experiment} is chosen - new experiment')
    proteins_predict_new_experiment_button.setVisible(True)
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
# end new


@magicgui(call_button='Create Segmentation')
def create_seggmentation():
    show_info("starting to create segmentations")
    segmentation.main(root_directory_path)
    show_info('done to create segmentations')
    return


@magicgui(call_button='New Experiment')
def new_exp():
    create_segmentation_button.setVisible(True)
    create_CSV_button.setVisible(True)

    upload_segmentation_button.setVisible(False)
    upload_csv_button.setVisible(False)
    patient_selection_button.setVisible(False)
    #upload_images_button.setVisible(False)
    choose_model_button.setVisible(False)
    ranking_model_button.setVisible(False)
    k_proteins_predict_button.setVisible(False)
    protein_selection_button.setVisible(False)
    change_std_button.setVisible(False)
    upload_csv_new_experiment_button.setVisible(False)

    new_exp_button.setVisible(False)
    old_exp_button.setVisible(True)
    return


@magicgui(call_button='Old Experiment')
def old_exp():
    upload_segmentation_button.setVisible(True)
    upload_csv_button.setVisible(True)

    create_segmentation_button.setVisible(False)
    create_CSV_button.setVisible(False)
    new_exp_button.setVisible(True)
    old_exp_button.setVisible(False)
    return


@magicgui(call_button='Upload Experiment')
def upload_exp():
    root = tk.Tk()
    root.withdraw()

    global root_directory_path
    try:
        # open the directory dialog box and allow the user to select a directory
        root_directory_path = filedialog.askdirectory(title="choose the root experiment directory")
        show_info(f'root experiment directory chosen successfully')
    except:
        show_info("add path to root experiment directory")
    upload_exp_button.setVisible(False)
    new_exp_button.setVisible(True)
    old_exp_button.setVisible(True)
    return


@magicgui
def my_text_box(my_text: str = "Enter text here") -> None:
    text_box = magicgui.widgets.LineEdit(text=my_text)
    return text_box


DEFAULT_CHOICES = []
DEFAULT_CHOICES_find_anomaly = []
DEFAULT_CHOICES_PATIENT = []


def _update_choices_on_file_change(widget):
    """Return available choices of proteins from csv headers"""
    filename = widget.filename.value
    choices = None
    choices_patient = None
    choices_find_anomaly = None
    global df_new_experiment
    widget.dropdown.choices = choices = ()
    widget.dropdown_patient.choices = ()
    widget.dropdown_find_anomaly.choices = ()
    if filename.is_file():
        if filename.suffix == ".csv":
            df_new_experiment = pd.read_csv(filename)
            choices = sorted(set(df.columns) - set(df_new_experiment.columns))
            choices_patient = df_new_experiment["SampleID"].drop_duplicates()
            choices_find_anomaly = np.setdiff1d(np.intersect1d(df.columns, df_new_experiment.columns),np.array(['SampleID', 'cellLabelInImage', 'cellSize']))
            print(choices_patient)
            print(choices)
            print(choices_find_anomaly)
        else:
            show_warning(f"File {filename} is not a .csv file.")
    if choices is not None:
        widget.dropdown.choices = choices
    if choices_patient is not None:
        print(type(widget.dropdown_patient.choices))
        print(widget.dropdown_patient.choices)
        widget.dropdown_patient.choices = choices_patient
        print(tuple(sorted(widget.dropdown_patient.choices)))
        widget.dropdown_patient.choices = tuple(sorted(widget.dropdown_patient.choices))
        print(widget.dropdown_patient.choices)
    if choices_find_anomaly is not None:
        widget.dropdown_find_anomaly.choices = choices_find_anomaly

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
        dropdown_find_anomaly=dict(
            widget_type="ComboBox", choices=DEFAULT_CHOICES_find_anomaly, label="Proteins to predict- Find anomaly",
        ),
        other_button=dict(widget_type="PushButton", text="Find Anomaly"),
        call_button="Predict Proteins",
        
    )

    def widget(viewer: napari.Viewer, dropdown_patient, dropdown,dropdown_find_anomaly,other_button, filename=Path.home(), foldername=Path.home()):
    #def widget(dropdown, dropdown_patient, filename=Path.home()):
        # Perform the prediction here
        proteins_list_to_predict = dropdown
        patient_number = dropdown_patient
        print(proteins_list_to_predict)
        print(patient_number)
        patient_number_new_experiment = patient_number
        #proteins_list = [protein.name for protein in choose_Proteins_New_Experiment]

        #print(protein_prediction_options_new_exeriment)
        list_of_proteins_to_train = np.setdiff1d(np.intersect1d(df.columns, df_new_experiment.columns),np.array(['SampleID', 'cellLabelInImage', 'cellSize']))
        # if (len(proteins_list) == 0):
        #     show_info("please select proteins - new Experiment")
        #     return
        # if df_new_experiment_normalized is None:
        #     show_info("upload csv first")
        #     return
        # if patient_number_new_experiment is None:
        #     show_info("choose patient number first (New Experiment)")
        #     return
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
        patient_cellLabel_image = find_cell_labeled_image(image_paths, patient_number_new_experiment)



        show_info("starting to predict proteins New Experiment")
        list_of_proteins_to_predict = proteins_list_to_predict  # ["CD45", "dsDNA", "Vimentin"]
        print(patient_number_new_experiment)
        print(proteins_list_to_predict)
        print(list_of_proteins_to_train)
        new_experiment.predict_k_proteins(viewer, df,df_new_experiment,patient_number_new_experiment, proteins_list_to_predict, list_of_proteins_to_train, patient_cellLabel_image)
        show_info('done predict proteins')



    @widget.filename.changed.connect
    def update_choices_on_file_change(event=None):
        _update_choices_on_file_change(widget)

    @widget.other_button.clicked.connect
    def do_something_with_other_button():
        # Perform other run button function here
        print("Clicked the other button")
        print(widget.dropdown_find_anomaly.value)

        patient_number_new_experiment = widget.dropdown_patient.value
        protein_to_predict = widget.dropdown_find_anomaly.value
        print(protein_to_predict)
        print(patient_number_new_experiment)

        list_of_proteins_to_train = np.setdiff1d(np.intersect1d(df.columns, df_new_experiment.columns),['SampleID', 'cellLabelInImage', 'cellSize', protein_to_predict])
        print(list_of_proteins_to_train)
        new_experiment.find_anomaly(viewer, df, df_new_experiment, patient_number_new_experiment,
                                          protein_to_predict, list_of_proteins_to_train, patient_cellLabel_image)


    return widget


def find_cell_labeled_image(paths, patient_number):
    for path in paths:
        file_name = os.path.basename(path)
        print(file_name)

        match_1 = re.search(r"p(\d+)_labeledcellData\.tiff", file_name)
        match_2 = re.search(r"p(\d+)_labeledcellData\.tif", file_name)
        if match_1:
            img_number_1 = int(match_1.group(1))
            if img_number_1 == patient_number:
                print(patient_number)
                img = Image.open(path)
                array_img = np.asarray(img)
                return array_img
        if match_2:
            img_number_2 = int(match_2.group(1))
            if img_number_2 == patient_number:
                print(patient_number)
                img = Image.open(path)
                array_img = np.asarray(img)
                return array_img
    return None








text_box_button = viewer.window.add_dock_widget(my_text_box, area='right')

# upload, new, old experiment buttons:
upload_exp_button = viewer.window.add_dock_widget(upload_exp, area='right')
new_exp_button = viewer.window.add_dock_widget(new_exp, area='right')
old_exp_button = viewer.window.add_dock_widget(old_exp, area='right')
# create buttons:
create_segmentation_button = viewer.window.add_dock_widget(create_seggmentation, area='right')
create_CSV_button = viewer.window.add_dock_widget(create_CSV, area='right')

# uplode buttons:
upload_segmentation_button = viewer.window.add_dock_widget(upload_segmentation, area='right')
upload_csv_button = viewer.window.add_dock_widget(upload_csv, area='right')
choose_model_button = viewer.window.add_dock_widget(choose_model, area='right')
ranking_model_button = viewer.window.add_dock_widget(rankingg_model, area='right')
patient_selection_button = viewer.window.add_dock_widget(patient_selection, area='right')
protein_selection_button = viewer.window.add_dock_widget(protein_selection, area='right')
change_std_button = viewer.window.add_dock_widget(widget_demo, area='right')
k_proteins_predict_button = viewer.window.add_dock_widget(proteins_predict, area='right')


upload_csv_new_experiment_button = viewer.window.add_dock_widget(upload_csv_new_experiment, area='right')
patient_selection_new_experiment_button = viewer.window.add_dock_widget(patient_selection_new_experiment, area='right')
proteins_predict_new_experiment_button = viewer.window.add_dock_widget(proteins_predict_new_experiment, area='right')
viewer.window.add_dock_widget(proteins_predict1(), name="Predict proteins new experiment",area='right' )

new_exp_button.setVisible(False)
old_exp_button.setVisible(False)
create_segmentation_button.setVisible(False)
create_CSV_button.setVisible(False)
upload_segmentation_button.setVisible(False)
upload_csv_button.setVisible(False)
patient_selection_button.setVisible(False)
choose_model_button.setVisible(False)
ranking_model_button.setVisible(False)
k_proteins_predict_button.setVisible(False)
protein_selection_button.setVisible(False)
change_std_button.setVisible(False)

#TO DO: change to set visble False
upload_csv_new_experiment_button.setVisible(True)
patient_selection_new_experiment_button.setVisible(True)
proteins_predict_new_experiment_button.setVisible(True)

def message():
    show_info('Welcome to Napari Plugin')


def main():
    napari


if __name__ == "__main__":
    main()
