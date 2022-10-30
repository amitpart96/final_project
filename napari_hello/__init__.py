import napari
import napari_hello
from napari_hello import run
from napari_hello import create_csv
from napari_hello import csv_test
from napari.utils.notifications import show_info
from skimage.io import imread


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
