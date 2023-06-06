from deepcell.applications import Mesmer
from skimage.io import imread
import numpy as np
from PIL import Image
import tkinter as tk
from tkinter import filedialog
import os


def save_img(seg, file_name):
    """
    saves a NumPy array representation of an image (seg) as a TIFF file with the specified file name.

    Parameters:
    - seg (NumPy array): A NumPy array representing the image to be saved.
    - file_name (str): The name of the file to be saved.

    Returns:
    - image_filename (str): The name of the saved image file.
    """
    new_im = Image.fromarray(seg)
    image_filename = f'{file_name}.tiff'
    new_im.save(image_filename)
    return image_filename


def main(root_directory_path):
    root = tk.Tk()
    root.withdraw()
    root_dir = root_directory_path

    # find the subfolders of the patients - each subbfolder is one patient that contains his proteins and a segmentation
    list_subfolders_with_paths = [f.path for f in os.scandir(root_dir) if f.is_dir()]
    list_patient = [f.name for f in os.scandir(root_dir) if f.is_dir()]
    subfiles_patient1 = list_subfolders_with_paths[0]

    # choose nuclear
    root = tk.Tk()
    root.withdraw()
    root_dir1 = filedialog.askopenfiles(initialdir=subfiles_patient1, title="choose images for for nuclear")
    filespath = [os.path.splitext(f.name)[0] for f in root_dir1]
    filesname_nuclear = [os.path.basename(filename) for filename in
                         filespath]  # get the base filename without the directory path

    # choose membrane
    root = tk.Tk()
    root.withdraw()
    root_dir1 = filedialog.askopenfiles(initialdir=subfiles_patient1, title="choose images for for membrane")
    filespath = [os.path.splitext(f.name)[0] for f in root_dir1]
    filesname_membrane = [os.path.basename(filename) for filename in
                          filespath]  # get the base filename without the directory path
    path_for_some_image =  str(root_dir) + "/" + list_patient[0] + "/" + filesname_nuclear[0] + ".tiff"
    some_image = np.asarray(Image.open(path_for_some_image))

    for patient in list_patient:
        path = str(root_dir) + "/" + patient
        sum_nuclear = np.zeros((np.size(some_image,0), np.size(some_image,1)), dtype=np.uint8)
        sum_membrane = np.zeros((np.size(some_image,0), np.size(some_image,1)), dtype=np.uint8)
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
        seg = labeled_image[0, :, :, 0]
        save_img(seg, path + "/SegmentationInterior")


if __name__ == "__main__":
    main(root_directory_path)


