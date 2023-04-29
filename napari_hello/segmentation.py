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


def main(root_directory_path):
    root = tk.Tk()
    root.withdraw()
    # root_dir = filedialog.askdirectory(title="Select the root folder containing all the patients folders")
    root_dir = root_directory_path
    print(root_dir)

    # find the subfolders of the patients - each sujbfolder is one patient that contains his proteins and a segmantation
    list_subfolders_with_paths = [f.path for f in os.scandir(root_dir) if f.is_dir()]
    print(list_subfolders_with_paths)
    list_patient = [f.name for f in os.scandir(root_dir) if f.is_dir()]
    subfiles_patient1 = list_subfolders_with_paths[0]

    # choose nuclear
    root = tk.Tk()
    root.withdraw()
    root_dir1 = filedialog.askopenfiles(initialdir=subfiles_patient1, title="choose images for for nuclear")
    filespath = [os.path.splitext(f.name)[0] for f in root_dir1]
    filesname_nuclear = [os.path.basename(filename) for filename in
                         filespath]  # get the base filename without the directory path
    print("choose nuclear proteins: ", filesname_nuclear)

    # choose membrane
    root = tk.Tk()
    root.withdraw()
    root_dir1 = filedialog.askopenfiles(initialdir=subfiles_patient1, title="choose images for for membrane")
    filespath = [os.path.splitext(f.name)[0] for f in root_dir1]
    filesname_membrane = [os.path.basename(filename) for filename in
                          filespath]  # get the base filename without the directory path
    print("choose membrane proteins: ", filesname_membrane)

    for patient in list_patient:
        path = root_dir + "/" + patient
        sum_nuclear = np.zeros((2048, 2048), dtype=np.uint8)
        sum_membrane = np.zeros((2048, 2048), dtype=np.uint8)
        for protein in filesname_nuclear:
            sum_nuclear += np.asarray(Image.open(path + "/" + protein + ".tiff"))
        for protein in filesname_membrane:
            sum_membrane += np.asarray(Image.open(path + "/" + protein + ".tiff"))

        # Combined together and expand to 4D
        im = np.stack((sum_nuclear, sum_membrane), axis=-1)
        im = np.expand_dims(im, 0)

        # create the application
        app = Mesmer()

        # create the lab
        labeled_image = app.predict(im)
        print(type(labeled_image))
        print(labeled_image.shape)
        labeled_image_reshape = labeled_image.reshape(2048, 2048)
        # labeled_image_reshape[labeled_image_reshape > 0] = 1
        seg = labeled_image[0, :, :, 0]
        print(seg.shape)
        save_img(seg, path + "/SegmentationInterior")


if __name__ == "__main__":
    main(root_directory_path)
