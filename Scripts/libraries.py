import ipyvolume as ipv
from ipywidgets import interact, interact_manual, widgets
import pylab as pl
from collections import Counter

# Biblioteca para estandarización de los datos
from numpy import asarray
from sklearn.preprocessing import MinMaxScaler

# Bibliotecas desarrolladas para el proyecto
# Juego del Caos y Dimensión Fractal
import os
import sys
#sys.path.insert(0, "../../DROGC-Bastida/Scripts/")

# Bilbioteca desarrollada para el presente proyecto
import SFI, DNA, fractal_dimension
from DNA import *
from SFI import *

#from fractal_dimension import *

"""
Para una graficación con diferente color de fondo
"""
#plt.style.use('dark_background')