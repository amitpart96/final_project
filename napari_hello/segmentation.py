from deepcell.applications import Mesmer
from skimage.io import imread
import numpy as np
from PIL import Image
import tkinter as tk
from tkinter import filedialog
import os

def save_img(seg, file_name):
    new_im = Image.fromarray(seg)
    # new_im.show()
    image_filename = f'{file_name}.tiff'
    print(image_filename)
    # save image using extension
    new_im.save(image_filename)
    return image_filename


def main():
    root = tk.Tk()
    root.withdraw()
    root_dir = filedialog.askdirectory()
    print(root_dir)

    # find the subfolders of the patients - each sujbfolder is one patient that contains his proteins and a segmantation
    list_subfolders_with_paths = [f.path for f in os.scandir(root_dir) if f.is_dir()]
    print(list_subfolders_with_paths)
    subfiles_patient1 = list_subfolders_with_paths[0]

    #choose nuclear
    root = tk.Tk()
    root.withdraw()
    root_dir1 = filedialog.askopenfiles(initialdir=subfiles_patient1,title="choose images for for nuclear")
    filespath = [os.path.splitext(f.name)[0] for f in root_dir1]
    filesname_nuclear = [os.path.basename(filename) for filename in filespath]  # get the base filename without the directory path
    print ("choose nuclear proteins: " , filesname_nuclear)

    #choose membrane
    root = tk.Tk()
    root.withdraw()
    root_dir1 = filedialog.askopenfiles(initialdir=subfiles_patient1, title="choose images for for membrane")
    filespath = [os.path.splitext(f.name)[0] for f in root_dir1]
    filesname_membrane = [os.path.basename(filename) for filename in
                         filespath]  # get the base filename without the directory path
    print("choose membrane proteins: ", filesname_membrane)

    all_proteins = [os.path.splitext(f.name)[0] for f in os.scandir(subfiles_patient1) if os.path.splitext(f.name)[1] == '.tiff']
    print(all_proteins)




    dsDNA = Image.open("C:/Users/amitp/Desktop/Point1/dsDNA.tiff")
    H3K27me3 = Image.open("C:/Users/amitp/Desktop/Point1/H3K27me3.tiff")
    H3K9ac = Image.open("C:/Users/amitp/Desktop/Point1/H3K9ac.tiff")
    dsDNA_Array = np.asarray(dsDNA)
    H3K27me3_Array = np.asarray(H3K27me3)
    H3K9ac_Array = np.asarray(H3K9ac)
    sum_nuclear = dsDNA_Array + H3K9ac_Array + H3K27me3_Array

    CD45 = Image.open("C:/Users/amitp/Desktop/Point1/CD45.tiff")
    HLADR = Image.open("C:/Users/amitp/Desktop/Point1/HLA-DR.tiff")
    Pan_keratin = Image.open("C:/Users/amitp/Desktop/Point1/Pan-keratin.tiff")
    beta_catenin = Image.open("C:/Users/amitp/Desktop/Point1/beta catenin.tiff")
    CD45_Array = np.asarray(CD45)
    HLADR_Array = np.asarray(HLADR)
    Pan_keratin_Array = np.asarray(Pan_keratin)
    beta_catenin_Array = np.asarray(beta_catenin)
    sum_membrane = CD45_Array + HLADR_Array + Pan_keratin_Array + beta_catenin_Array

    im1 = imread("C:/Users/amitp/Desktop/Point1/Vimentin.tiff")
    im2 = imread("C:/Users/amitp/Desktop/Point1/dsDNA.tiff")

    # Combined together and expand to 4D
    im = np.stack((sum_nuclear, sum_membrane), axis=-1)
    im = np.expand_dims(im, 0)

    # create the application
    app = Mesmer()

    # create the lab
    labeled_image = app.predict(im)
    print(type(labeled_image))
    print(labeled_image.shape)
    labeled_image_reshape = labeled_image.reshape(2048,2048)
    labeled_image_reshape[labeled_image_reshape > 0] = 1
    seg = labeled_image[0,:,:,0]
    print(seg.shape)
    save_img(labeled_image_reshape,"segmentation_amitAndLidor2")

if __name__ == "__main__":
    main()