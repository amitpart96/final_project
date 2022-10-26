import math
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
import csv
import plotly
import plotly.express as px
import plotly.graph_objects as go
from napari.utils.notifications import show_info
from skimage import data, filters, measure, morphology



def main():
    show_info('Hello, main!')
    oldIMG = Image.open('napari_hello/GtSegmentationInterior.tiff')
    #convert img to np array
    oldIMG = np.array(oldIMG)

    #convert img to binary img
    img = np.where(oldIMG ==0 , 0, 1)


    # Binary image, post-process the binary mask and compute labels
    threshold = filters.threshold_otsu(img)
    mask = img > threshold
    mask = morphology.remove_small_objects(mask, 100)
    mask = morphology.remove_small_holes(mask, 100)
    labels = measure.label(mask)

    fig = px.imshow(img, binary_string=True)
    fig.update_traces(hoverinfo='skip')  # hover is only for label info



    props = measure.regionprops(labels, img)
    properties = ['area', 'eccentricity', 'perimeter', 'intensity_mean']


    # For each label, add a filled scatter trace for its contour,
    # and display the properties of the label in the hover of this trace.
    for index in range(1, labels.max()):
        label_i = props[index].label
        contour = measure.find_contours(labels == label_i, 0.5)[0]
        y, x = contour.T
        hoverinfo = ''
        for prop_name in properties:
            hoverinfo += f'<b>{prop_name}: {getattr(props[index], prop_name):.2f}</b><br>'
        fig.add_trace(go.Scatter(
            x=x, y=y, name=label_i,
            mode='lines', fill='toself', showlegend=False,
            hovertemplate=hoverinfo, hoveron='points+fills'))


    plotly.io.show(fig)
    print("done show")




# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

