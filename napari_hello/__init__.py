import numpy as np

import napari_hello
from napari_hello import run
from napari_hello import create_csv
from napari_hello import csv_test
from tkinter import *
import napari
import napari_hello
from napari.types import ImageData, LayerDataTuple
from magicgui import magicgui
from enum import Enum
from napari.utils.notifications import show_info
from napari_hello import ranking_model


class Options(Enum):
    Au = 'Au'
    B7H3 = 'B7H3'
    Background = 'Background'
    Beta = 'Beta'
    catenin = 'catenin'
    C = 'C'
    Ca = 'Ca'
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
    CD163 = 'CD163'
    CD209 = 'CD209'
    CellTypes = 'CellTypes'
    CSF_1R = 'CSF-1R'
    dsDNA = 'dsDNA'
    EGFR = 'EGFR'
    Fe = 'Fe'
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
    Na = 'Na'
    OX40 = 'OX40'
    p = 'p'
    p53 ='p53'
    Pan_Keratin = 'Pan-Keratin'
    PD1 = 'PD1'
    PD_L1 = 'PD-L1'
    phospho_S6 = 'phospho-S6'
    Si = 'Si'
    SMA = 'SMA'
    Ta = 'Ta'
    Vimentin = 'Vimentin'


@magicgui(call_button='Create CSV')
def csv():
    create_csv.main()
    show_info('create csv')
    return

@magicgui(call_button='Ranking model')
def rankingg_model():
    ranking_model.main()
    show_info('ranking model')
    return

viewer = napari.Viewer()
viewer.window.add_dock_widget(csv, area='right')
viewer.window.add_dock_widget(rankingg_model, area='right')


@magicgui()
def protein_selection(protein_selection: Options) -> LayerDataTuple:
    # Do something with image and list of selected options
    print(protein_selection.value)
    show_info(protein_selection.value)
    return


viewer.window.add_dock_widget(protein_selection)






# 'Au','B7H3', 'Background', 'Beta catenin', 'C', 'Ca', 'CD3', 'CD4', 'CD8', 'CD11b', 'CD11c', 'CD16', 'CD20' ,
# 'CD31', 'CD45', 'CD45RO', 'CD56', 'CD63', 'CD68', 'CD68', 'CD163', 'CD209', 'CD209', 'CellTypes', 'CSF-1R',
# 'dsDNA','EGFR','Fe','FoxP3','H3K9ac', 'H3K27me3','HLA_class_1','HLA-DR','IDO','Keratin6','Keratin17','Ki67','Lag3',
# 'MPO','Na','OX40','p','p53','Pan-Keratin','PD1','PD-L1', 'phospho-S6','Si','SMA','Ta','Vimentin'


def main():
    napari


if __name__ == "__main__":
    main()


def show_hello_message():
    show_info('Hello, world!')


def csv():
    create_csv.main()
    show_info('Hello, run!')


def csv_test_func():
    csv_test.main()


# def ranking_model():
#     ranking_model.main()
