########### Generic
import numpy as np
import pandas as pd
from numpy import random

# Visualization libraries
import matplotlib.pyplot as plt
from matplotlib import cbook
from matplotlib.colors import LightSource

from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

########### Dynamical Systems
import pynamical
from pynamical import simulate, bifurcation_plot, save_fig, logistic_map, phase_diagram
import pandas as pd, numpy as np, IPython.display as display, matplotlib.pyplot as plt, matplotlib.cm as cm
from IPython.display import Markdown

# Ajuste de curva
from lmfit import Model

# Importando las bibliotecas para la interacción con las variables
#from __future__ import print_function
from ipywidgets import interact, interactive, fixed, interact_manual
import ipywidgets as widgets


########## Juego del Caos ##########

def resta_vectores_lista(A,B):
    return [i-j for i,j in zip(A,B)]
    
def multiply_lista_factor(factor,A):
    return [factor * i for i in A]
    
def suma_vectores_lista(A,B):
    return [i+j for i,j in zip(A,B)]

def rotacion_2d(theta, vector):
    r = np.array(( (np.cos(theta), -np.sin(theta)),
               (np.sin(theta),  np.cos(theta)) ))

    return r.dot(vector)


def regular_polygon(N,radio=1,anclaje=[0,0]):
    X_coordinate = []
    Y_coordinate = []
    Coordinate_Dictionary = {}    
    angle = 0

    angle_increment = 2 * np.pi / N
    for iteracion in range(N):
        X_coordinate.append( anclaje[0] + radio * np.cos(angle) )
        Y_coordinate.append( anclaje[1] + radio * np.sin(angle) )
        Coordinate_Dictionary[iteracion+1] = [X_coordinate[-1],Y_coordinate[-1]]
        angle += angle_increment

    return X_coordinate, Y_coordinate, Coordinate_Dictionary

def rotate_polygon(angle,polygon):
      return None
    
def norma_vector(vector):
    return np.linalg.norm(vector)


def chaos_game_polygon(N=10000,k=1/3,
              Polygon = {1:[0,0],
              2:[1,0],
              3:[1,1],
              4:[0,1]}):

    # Número de iteraciones
    #N = 100000
    iteracion = 0

    # Lista donde guardar iteraciones
    equis = []
    ye = []

    # Definimos el punto inicial
    Initial = Polygon[np.random.randint(1,len(Polygon)+1)]
    equis.append(Initial[0])
    ye.append(Initial[1])

    while iteracion < N:
        Punto_Local = np.random.randint(1,len(Polygon)+1)

        if Polygon[Punto_Local] == Initial:
            pass
        else:
            if type(k) == np.float:
                Local_Vector_Diff = multiply_lista_factor(k,
                                    resta_vectores_lista(Polygon[Punto_Local],Initial))

            else:
                k_ = random.choice(k)
                Local_Vector_Diff = multiply_lista_factor(k_,
                                    resta_vectores_lista(Polygon[Punto_Local],Initial))
                
            New_Position = suma_vectores_lista(Local_Vector_Diff, Polygon[Punto_Local])

            equis.append(New_Position[0])
            ye.append(New_Position[1])

            Initial = New_Position
        iteracion += 1
        
    return equis,ye



def chaos_game_polygon_(N=10000,k=[1/3],
              Polygon = {1:[0,0],
              2:[1,0],
              3:[1,1],
              4:[0,1]}):
    
    assert type(k) == list, "Los factores de escala deben estar dados por una lista, aunque sea un único elemento"

    #N = 100000
    iteracion = 0

    # Lista donde guardar iteraciones
    equis = []
    ye = []

    # Definimos el punto inicial
    Initial = Polygon[np.random.randint(1,len(Polygon)+1)]
    equis.append(Initial[0])
    ye.append(Initial[1])

    while iteracion < N:
        Punto_Local = np.random.randint(1,len(Polygon)+1)

        #if type(k) == np.float:
        #Local_Vector_Diff = multiply_lista_factor(k,resta_vectores_lista(Polygon[Punto_Local],Initial))

        #else:
        k_ = random.choice(k)
        Local_Vector_Diff = multiply_lista_factor(k_,resta_vectores_lista(Polygon[Punto_Local],Initial))
        New_Position = suma_vectores_lista(Initial, Local_Vector_Diff)

        equis.append(New_Position[0])
        ye.append(New_Position[1])

        Initial = New_Position
        iteracion += 1
        
    return equis,ye

def chaos_game_polygon_i(N=10000,
                         k=[1/3],
                         Polygon = {1:[0,0],2:[1,0],3:[1,1],4:[0,1]},
                         i = 1):
    
    assert type(k) == list, "Los factores de escala deben estar dados por una lista, aunque sea un único elemento"

    #N = 100000
    iteracion = 0

    # Lista donde guardar iteraciones
    equis = []
    ye = []

    # Definimos el punto inicial
    Initial = Polygon[i]
    equis.append(Initial[0])
    ye.append(Initial[1])

    while iteracion < N:
        Punto_Local = np.random.randint(1,len(Polygon)+1)

        k_ = random.choice(k)
        Local_Vector_Diff = multiply_lista_factor(k_,resta_vectores_lista(Polygon[Punto_Local],Initial))
        New_Position = suma_vectores_lista(Initial, Local_Vector_Diff)

        equis.append(New_Position[0])
        ye.append(New_Position[1])

        Initial = New_Position
        iteracion += 1

        i = Punto_Local
        
    return equis,ye

def chaos_game_polygon_i_no_rep(N=10000,
                                k=[1/3],
                                Polygon = {1:[0,0],2:[1,0],3:[1,1],4:[0,1]},
                                i = 1):
    
    assert type(k) == list, "Los factores de escala deben estar dados por una lista, aunque sea un único elemento"

    #N = 100000
    iteracion = 0

    # Lista donde guardar iteraciones
    equis = []
    ye = []

    # Definimos el punto inicial
    Initial = Polygon[i]
    equis.append(Initial[0])
    ye.append(Initial[1])

    while iteracion < N:
        Punto_Local = np.random.randint(1,len(Polygon)+1)

        if Punto_Local == i:
            pass
        #Local_Vector_Diff = multiply_lista_factor(k,resta_vectores_lista(Polygon[Punto_Local],Initial))

        else:
            k_ = random.choice(k)
            Local_Vector_Diff = multiply_lista_factor(k_,resta_vectores_lista(Polygon[Punto_Local],Initial))
            New_Position = suma_vectores_lista(Initial, Local_Vector_Diff)

            equis.append(New_Position[0])
            ye.append(New_Position[1])

            Initial = New_Position
            iteracion += 1
            
            i = Punto_Local
        
    return equis,ye


def Sier_few_iterations(random_state=124):
    P_N_triangle = 3
    Poli = regular_polygon(P_N_triangle)
    Legend = [2,3,1]
    
    Poligono_Rotado_Graficacion = [[],[]]
    Poli_rotado_ChaosGame = {i:None for i in list(Poli[2].keys())}

    # Rotamos el conjunto para alienarlo
    # Añadimos un contador para las etiquetas de la lista Legend
    indice = 0
    for vertice in range(P_N_triangle):
        Rotado_local = rotacion_2d(-np.pi/6,np.array((Poli[0][vertice],Poli[1][vertice])))
        Poligono_Rotado_Graficacion[0].append(Rotado_local[0])
        Poligono_Rotado_Graficacion[1].append(Rotado_local[1])

        Poli_rotado_ChaosGame[Legend[indice]] = [Rotado_local[0],Rotado_local[1]]
        indice += 1
    
    
    #plt.scatter(Poli[0], Poli[1])
    np.random.seed(random_state)

    fig=plt.figure(figsize=(9, 9), dpi= 80, facecolor='w', edgecolor='k')

    plt.scatter(Poligono_Rotado_Graficacion[0], Poligono_Rotado_Graficacion[1],s=70,c='k')

    for i, txt in enumerate(Legend):
        plt.annotate(txt, 
                     (Poligono_Rotado_Graficacion[0][i]+0.05, Poligono_Rotado_Graficacion[1][i]+0.075),
                     size=50)

    plt.annotate(r"$p_{0}$", 
                 (Poligono_Rotado_Graficacion[0][i]+0.05, Poligono_Rotado_Graficacion[1][i]-0.065),
                 color='red',
                 size = 50)


    # Generamos unas cuantas iteraciones del juego del caos Sierpinski
    POLI_FEW_REPS = chaos_game_polygon_i(N=5, k=[1/2], Polygon=Poli_rotado_ChaosGame, i=1)
    plt.scatter(POLI_FEW_REPS[0][1:],POLI_FEW_REPS[1][1:],s=20,c='r')

    LISTA_PUNTOS = []
    # Generamos las etiquetas que usaremos para describir el proceso
    for i in range(1,len(POLI_FEW_REPS[0][1:])+1):
        LISTA_PUNTOS.append(r"p_{" + r"{}".format(i) + "}")

    for i, txt in enumerate(LISTA_PUNTOS):
        plt.annotate(r"${}$".format(txt), #r"$\displaystyle {}$".format(txt), 
                     (POLI_FEW_REPS[0][1:][i]+0.05, POLI_FEW_REPS[1][1:][i]+0.05),
                     color='red',
                     size = 55)

    plt.axis('square')
    plt.axis('off')



def Sierpinski(Iteraciones=100000):
    # Generando los vértices
    Puntos_triangle = {1:[0,0],
              2:[1,0],
              3:[1/2,np.cos(np.pi / 3)]}

    Puntos_triangle_equis = [i[0] for i in list(Puntos_triangle.values())]
    Puntos_triangle_ye = [i[1] for i in list(Puntos_triangle.values())]

    # Fijamos el factor por el cual queremos multiplicar
    k = 1/2

    # Número de iteraciones
    N = Iteraciones
    iteracion = 0

    # Lista donde guardar iteraciones
    equis = []
    ye = []

    # Definimos el punto inicial
    #np.random.seed(500)
    Initial = Puntos_triangle[np.random.randint(1,len(list(Puntos_triangle))+1)]
    equis.append(Initial[0])
    ye.append(Initial[1])
    
    while iteracion < N:
        Punto_Local = np.random.randint(1,len(list(Puntos_triangle))+1)

        #if Puntos[Punto_Local] == Initial:
        #    pass
        #else:
        Local_Vector_Diff = multiply_lista_factor(k,
                            resta_vectores_lista(Puntos_triangle[Punto_Local],Initial))

        New_Position = suma_vectores_lista(Initial, Local_Vector_Diff)

        equis.append(New_Position[0])
        ye.append(New_Position[1])

        Initial = New_Position

        iteracion += 1

    plt.scatter(equis,ye,s=0.0002,c='k')
    plt.scatter(Puntos_triangle_equis,Puntos_triangle_ye,c='k', s=5)
    #plt.axis('square')
    plt.axis('off')
    plt.show()

    
def Sierpinski_coordinates(Iteraciones=100000):
    # Generando los vértices
    Puntos_triangle = {1:[0,0],
              2:[1,0],
              3:[1/2,np.cos(np.pi / 3)]}

    Puntos_triangle_equis = [i[0] for i in list(Puntos_triangle.values())]
    Puntos_triangle_ye = [i[1] for i in list(Puntos_triangle.values())]

    # Fijamos el factor por el cual queremos multiplicar
    k = 1/2

    # Número de iteraciones
    N = Iteraciones
    iteracion = 0

    # Lista donde guardar iteraciones
    equis = []
    ye = []

    # Definimos el punto inicial
    #np.random.seed(500)
    Initial = Puntos_triangle[np.random.randint(1,len(list(Puntos_triangle))+1)]
    equis.append(Initial[0])
    ye.append(Initial[1])
    
    while iteracion < N:
        Punto_Local = np.random.randint(1,len(list(Puntos_triangle))+1)

        #if Puntos[Punto_Local] == Initial:
        #    pass
        #else:
        Local_Vector_Diff = multiply_lista_factor(k,
                            resta_vectores_lista(Puntos_triangle[Punto_Local],Initial))

        New_Position = suma_vectores_lista(Initial, Local_Vector_Diff)

        equis.append(New_Position[0])
        ye.append(New_Position[1])

        Initial = New_Position

        iteracion += 1

    #plt.scatter(equis,ye,s=0.0002,c='k')
    #plt.scatter(Puntos_triangle_equis,Puntos_triangle_ye,c='k', s=5)
    #plt.axis('square')
    #plt.axis('off')
    #plt.show()
    return equis, ye

#def henon(iterations,a=1.4,b=0.3,x0=0.5):
#    """
#    Calcula las coordenadas del mapa de henón
#    """
#    ye = [b*x0]
#    equis = [1 - a*(x0**2) + ye[-1]]
#   
#    for i in range(iterations-1):
#        ye.append(b*equis[-1])
#        equis.append(1 - a*(equis[-1]**2) + ye[-1])
#    return equis,ye


def henon_attractor(x, y, a, b):
    '''Computes the next step in the Henon 
    map for arguments x, y with kwargs a and
    b as constants.
    '''
    x_next = 1 - a * x ** 2 + y
    y_next = b * x
    return x_next, y_next

def henon(iteration, a=1.4, b=0.3):
    X = np.zeros(iteration + 1)
    Y = np.zeros(iteration + 1)

    # starting point
    X[0], Y[0] = 0, 0

    # add points to array
    for i in range(iteration):
        x_next, y_next = henon_attractor(X[i], Y[i], a, b)
        X[i+1] = x_next
        Y[i+1] = y_next
        
    return X,Y

    
def Chaos_Game(N=3,Iteraciones=10000,k1=1/2, k2=1/2):
    P_N = N
    Poli = regular_polygon(P_N)

    ITERACIONES = Iteraciones#100000
    FRACCIONES_RECORRIDAS = [k1,k2]
    Poli_Chaos = chaos_game_polygon_i(N=ITERACIONES, k=FRACCIONES_RECORRIDAS, Polygon=Poli[2])
    
    # Graficación
    fig=plt.figure(figsize=(9, 9), dpi= 80, facecolor='w', edgecolor='k')
    plt.scatter(Poli_Chaos[0],
                Poli_Chaos[1],s=0.002,c='k')

    plt.axis('square')
    plt.axis('off')
    plt.show()

    
################# Sistemas Dinámicos

## Variables
title_font = pynamical.get_title_font()
label_font = pynamical.get_label_font()

shapes = ['^','d','x','s','1','P',"o"]


# Visualization
def get_colors(cmap, n, start=0., stop=1., alpha=1., reverse=False):
    '''return n-length list of rgba colors from the passed colormap name and alpha,
       limit extent by start/stop values and reverse list order if flag is true'''
    colors = [cm.get_cmap(cmap)(x) for x in np.linspace(start, stop, n)]
    colors = [(r, g, b, alpha) for r, g, b, _ in colors]
    return list(reversed(colors)) if reverse else colors




####################### Para el ajuste de curva #############################
def Germi_Logaritmica_lmfit(x,K,b,N,a):
    return np.log(K) + b * np.log(N+1-x) - a * np.log(x)

def Germi_Beta_lmfit(x,K,b,a,N):
    return K * (((N + 1 - x) ** b) / (x ** a))
    #return K * ((N + 1 - x) ** b)



################### Miscelaneo ####################
def line(equis,a,b):
    return b + a * x

