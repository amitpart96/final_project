import time
import tkinter as tk
from idlelib.tooltip import _tooltip
from tkinter import filedialog as fd, filedialog

import PIL.Image
import pandas as pd
from PIL.Image import Image
from dask.array.chunk import view
from imageio import imread

import napari_hello
from napari_hello import create_csv
import napari
import napari_hello
import os
from magicgui import magicgui
from enum import Enum
from napari.utils.notifications import show_info
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
    Beta_catenin ='Beta catenin'
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

@magicgui(call_button='Create CSV')
def create_CSV():
    create_csv.main()
    show_info('created csv successfully')
    patient_selection_button.setVisible(True)
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
    find_anomaly.main(viewer, df, patient_number, protein)
    show_info('done find anomaly')
    return


@magicgui(call_button='Select Protein')
def protein_selection(protein_selection: Options_Proteins):
    # Do something with image and list of selected options
    global protein
    protein = protein_selection.value
    show_info(f'{protein} is chosen')
    find_anomaly_button.setVisible(True)
    return


@magicgui(call_button='Select Patient')
def patient_selection(patient_selection: Options_Patients):
    # Do something with image and list of selected options
    global patient_number
    patient_number = int(patient_selection.value)
    show_info(f'patient {patient_number} is chosen')
    ranking_model_button.setVisible(True)
    protein_selection_button.setVisible(True)
    return

@magicgui(call_button='Upload Images')
def upload_images():
    root = tk.Tk()
    root.withdraw()
    list_img = fd.askopenfilenames()  # choose celldata of the patient
    for img in list_img:
        napari_image = imread(img)  # Reads an image from file
        img_name=os.path.basename(img)
        if patient_number is not None:
            img_name="Patient" + str(patient_number) + " " + img_name
        viewer.add_image(napari_image, name=img_name)  # Adds the image to the viewer and give the image layer a name
    show_info(f'images uploaded successfully')
    return

upload_segmentation_button = viewer.window.add_dock_widget(upload_segmentation, area='right')
create_CSV_button = viewer.window.add_dock_widget(create_CSV, area='right')
upload_csv_button = viewer.window.add_dock_widget(upload_csv, area='right')
patient_selection_button = viewer.window.add_dock_widget(patient_selection, area='right')
ranking_model_button = viewer.window.add_dock_widget(rankingg_model, area='right')
protein_selection_button = viewer.window.add_dock_widget(protein_selection, area='right')
find_anomaly_button = viewer.window.add_dock_widget(findd_anomaly, area='right')
upload_images_button = viewer.window.add_dock_widget(upload_images, area='right')

patient_selection_button.setVisible(False)
upload_images_button.setVisible(False)
ranking_model_button.setVisible(False)
protein_selection_button.setVisible(False)
find_anomaly_button.setVisible(False)

def message():
    show_info('Welcome to Napari Plugin')


def main():
    napari


if __name__ == "__main__":
    main()
