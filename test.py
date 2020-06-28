from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
import pandas

data_frame = pandas.read_excel("results_bhe_300m.xlsx", sheet_name="results_bhe_300m")

x = data_frame["H_clay"]
y = data_frame["H_soil"]
z = data_frame["Q_extraction"]

fig = pyplot.figure()

ax = Axes3D(fig)

ax.scatter(x, y, z, "r.")

pyplot.show()
