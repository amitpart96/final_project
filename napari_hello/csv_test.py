import csv
import pandas as pd
from skimage import data, filters, measure, morphology
import tkinter as tk
from tkinter import filedialog

def csv():
    root = tk.Tk()
    root.withdraw()
    rootdir  = filedialog.askdirectory()
    print("1")
    print(rootdir)
    print("2")

    all_rows = []
    col_names_csv = ['patient Number', 'cell index', 'cell size', 'Au','B7H3', 'Background', 'Beta catenin', 'C', 'Ca', 'CD3', 'CD4', 'CD8', 'CD11b', 'CD11c', 'CD16', 'CD20' , 'CD31', 'CD45', 'CD45RO', 'CD56', 'CD63', 'CD68', 'CD68', 'CD163', 'CD209', 'CD209', 'CellTypes', 'CSF-1R', 'dsDNA','EGFR','Fe','FoxP3','H3K9ac','H3K27me3','HLA_class_1','HLA-DR','IDO','Keratin6','Keratin17','Ki67','Lag3','MPO','Na','OX40','p','p53','Pan-Keratin','PD1','PD-L1','phospho-S6','Si','SMA','Ta','Vimentin' ]
    all_rows.append(col_names_csv)
    df = pd.DataFrame(all_rows)
    df.to_csv('empty_csv.csv')


def main():
    csv()

