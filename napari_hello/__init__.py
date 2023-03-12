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


@magicgui(call_button='Upload Csv')
def upload_csv():
    root = tk.Tk()
    root.withdraw()
    try:
        global df
        filename = fd.askopenfilename()
        print(filename)
        df = pd.read_csv(filename)
        show_info(f'csv uploaded successfully')
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
    show_info('created csv successfully')
    create_CSV_button.setVisible(True)
    ranking_model_button.setVisible(True)
    patient_selection_button.setVisible(True)


def exepc_in_create_csv():
    create_CSV_button.setVisible(True)


@magicgui(call_button='Create CSV')
def create_CSV():
    # create_CSV_button.setVisible(False)
    show_info("Processing CSV creation")
    thread = threading.Thread(target=create_csv.main)
    thread.start()
    return


@magicgui(call_button='Rank Proteins')
def rankingg_model():
    if df is None:
        show_info("upload csv first")
        return
    if patient_number is None:
        show_info("choose patient number first")
        return
    show_info("starting to rank proteins")
    time.sleep(3)
    ranking_model.main(viewer, df, patient_number)
    show_info('done rank proteins')
    return


@magicgui(call_button='Find Anomaly')
def findd_anomaly():
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


@magicgui(call_button='Select Protein')
def protein_selection(protein_selection: Options_Proteins):
    # Do something with image and list of selected options
    global protein
    protein = protein_selection.value
    show_info(f'{protein} is chosen')
    find_anomaly_button.setVisible(True)
    change_std_button.setVisible(True)
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


# widget_demo.show()
upload_segmentation_button = viewer.window.add_dock_widget(upload_segmentation, area='right')
create_CSV_button = viewer.window.add_dock_widget(create_CSV, area='right')
upload_csv_button = viewer.window.add_dock_widget(upload_csv, area='right')
ranking_model_button = viewer.window.add_dock_widget(rankingg_model, area='right')
patient_selection_button = viewer.window.add_dock_widget(patient_selection, area='right')
protein_selection_button = viewer.window.add_dock_widget(protein_selection, area='right')
find_anomaly_button = viewer.window.add_dock_widget(findd_anomaly, area='right')
upload_images_button = viewer.window.add_dock_widget(upload_images, area='right')
change_std_button = viewer.window.add_dock_widget(widget_demo, area='right')
patient_selection_button.setVisible(False)
upload_images_button.setVisible(False)
ranking_model_button.setVisible(False)
protein_selection_button.setVisible(False)
find_anomaly_button.setVisible(False)
change_std_button.setVisible(False)


def message():
    show_info('Welcome to Napari Plugin')


def main():
    napari


if __name__ == "__main__":
    main()
