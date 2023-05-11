import shapefile as shp  # Requires the pyshp package
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection


def plot_island():
    sf = shp.Reader("/Users/pselvaraj/Github/testing/emodpy-vector_genetics/comoros_exploration/data/com_adm_cosep_ocha_20191205_shp/com_admbnda_adm1_cosep_ocha_20191205.shp")

    fig = plt.figure(figsize=(5, 7.5))
    ax = fig.add_subplot(111)
    for shape in sf.shapeRecords():
        if shape.record.ADM1_EN == 'Grande Comore (Ngazidja)':
            x = [i[0] for i in shape.shape.points[:]]
            y = [i[1] for i in shape.shape.points[:]]
            ax.plot(x, y, color='k', lw=0.5)

    plt.xticks([])
    plt.yticks([])
    plt.show()


if __name__ == '__main__':
    plot_island()