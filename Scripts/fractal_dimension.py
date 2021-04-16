import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from collections import Counter

def line(equis,a,b):
    return [b + a * x for x in equis]

def find_min_distance_list(LIST):
    A = sorted(LIST)
    RANGO = max(LIST) - min(LIST)
    B = [A[i] - A[i-1] for i in range(1,len(LIST))]
    B_min = min(B)
    Min_Index = next((i for i, x in enumerate(sorted(B)) if x), None)
    
    return B[Min_Index]

def range_finder(LISTA):
	return max(LISTA) - min(LISTA)


def log_2_scales(rango,minimo):
    criterio = False
    scales = []
    Initial = 1
    while not criterio:
        value_append = rango * (2**(-Initial))
        scales.append(value_append)
        Initial += 1
        if value_append < minimo:
            criterio = True
    
    return scales


def fractal_dimension(coord_list):
    MINIMUM_SCALE = max(find_min_distance_list(coord_list[0]),
                        find_min_distance_list(coord_list[1]))
    
    # Asumimos que la escala de una lista es suficiente porque la escala en ambos ejes es la misma
    scales = log_2_scales(range_finder(coord_list[0]),MINIMUM_SCALE)
    
    counter = []
    for scale in scales:
        counter.append(len([i for i in rec if i[0]<scale and i[1]<scale]))

    # Quitando los ceros 
    indice = next((i for i, x in enumerate(sorted(counter)) if x), None)
    if indice != 0:
        counter = counter[0:-indice]
        scales = scales[0:-indice]
    else:
        pass
    
    log_counter = [np.log(counter)]
    log_scales = np.log(1/np.array(scales))
    
    coeffs= np.polyfit(log_scales,log_counter[0], 1)
    print(f"coefficient (fractal dimension) = {np.abs(coeffs[0])}")
    
    plt.scatter(log_scales,log_counter[0])


############## Computation of the dimension through definition ###############

def bin_like(lista,escala):
	"""
	Partition of an interval given the scaling factor
	"""
	sorted_list = sorted(lista)
	next_ = escala
	result = [min(lista)]
	while result[-1] < max(lista):
	    result.append(result[-1] + next_)
	    next_ += escala
	    
	return result


def bin_like_basetwo(lista,N):
	"""
	Similar to bin with power 2 partitions
	"""
	scale = 2**(-N)
	lista = sorted(lista)
	Initial = [lista[0]]
	while Initial[-1] < max(lista):
	    Initial.append(Initial[-1]+scale)
	return Initial,scale



def Frac_dimension_computation(equis,ye):
	Intera = 2
	iteration = []
	EQUIS = []
	YE = []
	counter_iter = 1
	while bin_like_basetwo([min(equis),max(equis)],Intera)[1] > 1*find_min_distance_list(equis) and counter_iter < 8:
	    iteration.append(bin_like_basetwo([0,1],Intera)[0])
	    
	    # Histograma
	    bi = bin_like_basetwo([min(equis),max(equis)],Intera)[0]
	    H = np.histogram2d(equis,ye, bins=((bi,bi)))
	    YE.append(np.log(np.sum(H[0]>0)))
	    EQUIS.append(np.log(1/bin_like_basetwo([0,1],Intera)[1]))
	    
	    #print(bin_like_basetwo([0,1],Intera)[1],find_min_distance_list(equis))
	    #print(np.sum(H[0]>0))
	    
	    Intera += 1
	    counter_iter +=1


	coeffs= np.polyfit(EQUIS,YE, 1)
	print(f"coefficient (fractal dimension) = {np.abs(coeffs[0])}")

	plt.scatter(EQUIS,YE)
	plt.xlabel("log(1 / $\epsilon$)")
	plt.ylabel("log(N ($\epsilon$))")
	plt.show()
	plt.clf()

	return EQUIS, YE, coeffs[0]

def Frac_dimension_computation_name(equis,ye,name):
	Intera = 2
	iteration = []
	EQUIS = []
	YE = []
	counter_iter = 1
	while bin_like_basetwo([min(equis),max(equis)],Intera)[1] > 1*find_min_distance_list(equis) and counter_iter < 8:
	    iteration.append(bin_like_basetwo([0,1],Intera)[0])
	    
	    # Histograma
	    bi = bin_like_basetwo([min(equis),max(equis)],Intera)[0]
	    H = np.histogram2d(equis,ye, bins=((bi,bi)))
	    YE.append(np.log(np.sum(H[0]>0)))
	    EQUIS.append(np.log(1/bin_like_basetwo([0,1],Intera)[1]))
	    
	    #print(bin_like_basetwo([0,1],Intera)[1],find_min_distance_list(equis))
	    #print(np.sum(H[0]>0))
	    
	    Intera += 1
	    counter_iter +=1


	coeffs= np.polyfit(EQUIS,YE, 1)
	print(f"coefficient (fractal dimension) = {np.abs(coeffs[0])}")

	plt.scatter(EQUIS,YE)
	plt.xlabel("log(1 / $\epsilon$)")
	plt.ylabel("log(N ($\epsilon$))")
	plt.title(name)
	plt.savefig("{}".format(name),dpi=150)
	plt.show()
	plt.clf()

	return EQUIS, YE, coeffs[0]



def half_way(p1,p2):
    return ((p1[0]+p2[0])/2,(p1[1]+p2[1])/2)

def chaos_game_chromosome(string,coordinates_dic):
    resulting_coordinates = []
    INITIAL = coordinates_dic[string[0]]
    for i in string[1:]:
        resulting_coordinates.append(half_way(INITIAL,coordinates_dic[i]))
        INITIAL = resulting_coordinates[-1]
    return resulting_coordinates

def error(a,b):
	return abs(1 - (a/b))


#def random_nuc(N):
#	"""
#	Returns random nucleotide-like string
#	"""
#	d = {1:'A',2:'T',3:'C',4:'G'}
#	string = str()
#	for i in range(0,N):
#		string += d[np.random.randint(1,5)]

	return string


def cumulative_list(LIST):
	assert LIST == sorted(LIST), "Lista debe estar ordenada"
	cumulative_list = [[0,LIST[0]]]
	for i in range(1,len(LIST)):
		cumulative_list.append([1+cumulative_list[-1][1],cumulative_list[-1][1]+LIST[i]])
        
	return cumulative_list


def interval_picker(value,lista):
	assert value < lista[-1][1], "valor"
	crit = False
	initial = 0
	while not crit:
		if lista[initial][0] < value < lista[initial][1]:
			crit = True
		else:
			initial += 1
            
	return initial


def random_nuc(N,bias=False):
	"""
	Returns a biased random nucleotide-like string
	"""
	if bias == False:
		d = {1:'A',2:'T',3:'C',4:'G'}
		string = str()
		for i in range(0,N):
			string += d[np.random.randint(1,5)]

		return string
	else:
		sort_list = cumulative_list(list(sort_dict_values(invert_dic(bias),rev=False,by='keys').keys()))
		Rand = np.random.randint(1,sort_list[-1][1])
		string = str()
		for i in range(0,N):
			string += list(sort_dict_values(invert_dic(bias),rev=False,by='keys').values())[interval_picker(Rand,sort_list)]

		return string
    
def string_dist(string):
	"""
	Verifica la distribución que sigue una cadena de caracteres (string)
	"""
	string = list(string)
	return dict(Counter(string))



def interval_picker(value,lista):
	assert value < lista[-1][1], "valor"
	crit = False
	initial = 0
	while not crit:
		if lista[initial][0] < value < lista[initial][1]:
			crit = True
		else:
			initial += 1
	return initial


def invert_dic(dic):
	return dict(zip(dic.values(),dic.keys()))