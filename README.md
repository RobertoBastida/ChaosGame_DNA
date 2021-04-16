<center><h1>Dimensión fractal en diferentes conjuntos</h1></center>

El presente notebook tiene por objetivo ilustrar la ejecución de las funciones desarrolladas para el cálculo de la dimensión fractal para subconjunto de $\mathbb{R}^{2}$.


```python
import sys
sys.path.insert(0, "Scripts")

from libraries import *
from fractal_dimension import *
```

## Comparación entre conteo de cajas y dimensión topológica

Para verificar que las funciones se comporten correctamente; las probamos primero con dos conjuntos arquetípicos: una línea recta y un cuadrado. 

Consideramos válidas las funciones para el cálculo de la dimensión a partir del algoritmo de Minkowski-Bouling si es que éstas arrojaran las dimensiones iguales a 1.0 y 2.0 para los conjuntos antes mencionados, respectimante.


```python
length = 100000
```

## Recta

Generación de una recta con 1,000,000 puntos.


```python
line_x = np.linspace(0,1,length)
line_y = [0.5 for x in line_x]
```


```python
plt.scatter(line_x,line_y,s=0.00005)
plt.axis('equal')
plt.axis('off')
plt.savefig("ilustraciones/plot/Otros/{}".format("recta"),dpi=150)
```


![png](output_6_0.png)



```python
Line_dimension = Frac_dimension_computation_name(line_x,line_y,"ilustraciones/dim/Otros/Recta")
```

    coefficient (fractal dimension) = 1.0



![png](output_7_1.png)



    <Figure size 432x288 with 0 Axes>


## Cuadrado

Generación de puntos aleatorios dentro de un cuadrado


```python
# La generación de los puntos aleatorios
rec_coordinates = np.random.rand(length,2)

# Seleccionamos las x,y como listas
rec_x = rec_coordinates[:,0]
rec_y = rec_coordinates[:,1]
```


```python
# Graficamos para visualizar el conjunto
plt.scatter(rec_x,rec_y,s=0.00005)
plt.axis('equal')
plt.axis('off')
plt.savefig("ilustraciones/plot/Otros/{}".format("cuadrado"),dpi=150)
```


![png](output_10_0.png)



```python
Square_dimension = Frac_dimension_computation_name(rec_x,rec_y,"ilustraciones/dim/Otros/Cuadrado")
```

    coefficient (fractal dimension) = 1.9615362805847147



![png](output_11_1.png)



    <Figure size 432x288 with 0 Axes>


# Juego del caos

Probamos ahora el correcto funcionamiento de las librerías desarrolladas específicamente para medir la dimensión fractal de conjuntos generados a partir del juego del caos en dos diferentes escenarios: El triángulo de Sierpinski


```python
sier_x,sier_y = SFI.Sierpinski_coordinates(Iteraciones=length)
```


```python
plt.scatter(sier_x,sier_y,s=0.00005)
#plt.axis('equal')
plt.axis('off')
plt.savefig("ilustraciones/plot/Otros/{}".format("sierpinski"),dpi=150)
```


![png](output_14_0.png)



```python
Sierpinski_dimension = Frac_dimension_computation_name(sier_x,sier_y,"ilustraciones/dim/Otros/Sierpinski")
```

    coefficient (fractal dimension) = 1.6546925948509885



![png](output_15_1.png)



    <Figure size 432x288 with 0 Axes>



```python
line_r_x = sorted(sier_x)
line_r_y = [0.5 + 1*(np.random.rand()) for x in line_r_x]
```


```python
plt.scatter(line_r_x,line_r_y,s=0.00005)
plt.axis('equal')
plt.axis('off')
plt.savefig("ilustraciones/plot/Otros/{}".format("barcode"),dpi=150)
```


![png](output_17_0.png)



```python
Barcode_dimension = Frac_dimension_computation_name(line_r_x,line_r_y,"ilustraciones/dim/Otros/Código de barras")
```

    coefficient (fractal dimension) = 1.8581990300631681



![png](output_18_1.png)



    <Figure size 432x288 with 0 Axes>


## Genoma humano

Se emula el juego del caos con el genoma humano a basándonos en las ideas de <u>artículo</u>


```python
Dic_2 = {'A':(0,0),'C':(1,0),'G':(1,1),'T':(0,1)}
```


```python
Chr = "chr_19_.txt"
Str_Chr = load_txt_file("data/",Chr)

# trimming the string to a length of maximum 1,000,000 nucleotides
CG_Chr = chaos_game_chromosome(Str_Chr[0:min(len(Str_Chr),length)],Dic_2)

# grabbing the coordinates
CG_x = [i[0] for i in CG_Chr]
CG_y = [i[1] for i in CG_Chr]

# Plotting 
plt.scatter(CG_x,CG_y,s=0.0000125)
plt.axis('equal')
plt.axis('off')
plt.savefig("ilustraciones/plot/ch19/{}.png".format(Chr.split(".")[0]),dpi=200)
plt.show()
plt.clf()

Chr_19_Fd = Frac_dimension_computation_name(CG_x,CG_y,"ilustraciones/dim/ch19/{}".format(Chr.split(".")[0]))
```


![png](output_21_0.png)


    coefficient (fractal dimension) = 1.8975485424091463



![png](output_21_2.png)



    <Figure size 432x288 with 0 Axes>


## Random


```python
Str_Chr = random_nuc(length)

# trimming the string to a length of maximum 1,000,000 nucleotides
CG_Chr = chaos_game_chromosome(Str_Chr[0:min(len(Str_Chr),length)],Dic_2)

# grabbing the coordinates
CG_x = [i[0] for i in CG_Chr]
CG_y = [i[1] for i in CG_Chr]

plt.scatter(CG_x,CG_y,s=0.0000125)
plt.axis('equal')
plt.axis('off')
plt.savefig("ilustraciones/plot/random/{}.png".format("random"),dpi=200)
plt.show()
plt.clf()

Chr_19_Fd = Frac_dimension_computation_name(CG_x,CG_y,"ilustraciones/dim/random/{}".format("random"))
```


![png](output_23_0.png)


    coefficient (fractal dimension) = 1.9622509118831655



![png](output_23_2.png)



    <Figure size 432x288 with 0 Axes>

