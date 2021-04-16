import numpy as np
import scipy.stats
import pynamical
import matplotlib.pyplot as plt
import pandas as pd
from DNA import sort_dict_values

####### Variables #######

####### Visualization
shapes = ['^','d','x','s','1','P',"o"]
title_font = pynamical.get_title_font()
label_font = pynamical.get_label_font()

x_values = np.arange(-5, 5, 0.15)
x_scatter = np.linspace(-5,5,25)

####### Normal Distribution

def normal_plot_compare(mu,sigma,mu_1,sigma_1):
    global x_values
    y_values = scipy.stats.norm(mu, sigma)
    y_values_1 = scipy.stats.norm(mu_1, sigma_1)
    
    PARÁMETROS =  [[mu,sigma],[mu_1,sigma_1]]
    
    for param in PARÁMETROS:
        LOCAL_VALUES = scipy.stats.norm(param[0], param[1])
        plt.plot(x_values,LOCAL_VALUES.pdf(x_values),linewidth=2)
        plt.scatter(x_scatter,LOCAL_VALUES.pdf(x_scatter),marker=shapes[PARÁMETROS.index(param)],
                   label=r'$\mu$ = {}  $\sigma$ = {}'.format(str(round(param[0],2)),str(round(param[1],2))))
    
    plt.legend(title="Parámetros", loc='best', bbox_to_anchor=(1, 0.525))
    plt.xlabel('Equis', fontproperties=label_font)
    plt.ylabel('Probabilidad', fontproperties=label_font)
    plt.title('Normal', fontproperties=label_font)
    
    
####### Ley de Zipf #######
r = [i for i in range(1,21)]
Alphas = [0.25,0.50,0.75]

def Zipf_unsorted(r,alpha):
    #global r
    if type(r) == list:
        return [i**(-alpha) for i in r]
    elif type(r) in [int]:
        return r**(-alpha)
    
def Zipf(K,r,alpha):
    #global r
    return (K) / (r**alpha)

def Normalize_Zipf(r,alpha):
    #global r
    Anormal_list = Zipf_unsorted(r,alpha)
    #print(Anormal_list)
    SUM = sum(Anormal_list)
    K = SUM ** (-1)
    
    return [[Zipf(K,i,alpha) for i in r],r,alpha,K]

def Zipf_plot_compare(alpha,alpha_1):
    #global r
    PARÁMETROS =  [alpha,alpha_1]
    
    for param in PARÁMETROS:
        LOCAL_VALUES = Normalize_Zipf(r,param)
        plt.plot(r,LOCAL_VALUES[0],linewidth=2)
        plt.scatter(r,LOCAL_VALUES[0],marker=shapes[PARÁMETROS.index(param)],label=str(round(param,2)))
    
    plt.legend(title=r'$\alpha$', loc='best', bbox_to_anchor=(1, 0.525))
    plt.xlabel('Rango', fontproperties=label_font)
    plt.ylabel('Probabilidad', fontproperties=label_font)
    plt.title('Zipf', fontproperties=label_font)
    
    
    
####### MANDELBROT #######
#r = [i for i in range(1,21)]
Epsilon = [0.25,0.50,0.75]
Rho = [0.1,0.2,0.3]
r = [i for i in range(1,21)]

def Maldelbrot(r,epsilon,rho):
    N = max(r)
    return [((N + rho) / (i + N)) ** (1+epsilon) for i in r]

def Normalize_Maldelbrot(r,epsilon,rho):
    Anormal_list = Maldelbrot(r,epsilon,rho)
    #print(Anormal_list)
    SUM = sum(Anormal_list)
    K = SUM ** (-1)
    
    return [[K * i for i in Anormal_list],epsilon,rho]

def Maldelbrot_plot_compare(e_1,rho_1,e_2,rho_2):
    #global r
    PARÁMETROS =  [[e_1,rho_1],[e_2,rho_2]]
    
    for param in PARÁMETROS:
        LOCAL_VALUES = Normalize_Maldelbrot(r,param[0],param[1])
        plt.plot(r,LOCAL_VALUES[0],linewidth=2)
        plt.scatter(r,LOCAL_VALUES[0],marker=shapes[PARÁMETROS.index(param)],
                    label=r'$\epsilon$ = {}  $\rho$ = {}'.format(str(round(param[0],2)),str(round(param[1],2))))
    
    plt.legend(title="Parámetros", loc='best', bbox_to_anchor=(1, 0.525))
    plt.xlabel('Rango', fontproperties=label_font)
    plt.ylabel('Probabilidad', fontproperties=label_font)
    plt.title('Maldelbrot', fontproperties=label_font)
    
    
####### LAVALETTE #######
r = [i for i in range(1,21)]
Alphas = [0.25,0.50,0.75]

def Lavalette(r,b):
    N = max(r)
    return [((N + 1 - i) / (i)) ** (b) for i in r]

def Normalize_Lavalette(r,b):
    Anormal_list = Lavalette(r, b)
    #print(Anormal_list)
    SUM = sum(Anormal_list)
    K = SUM ** (-1)
    
    return [[K * i for i in Anormal_list],b]

def Lavalette_plot_compare(b,b_1):
    #global r
    PARÁMETROS =  [b,b_1]
    
    for param in PARÁMETROS:
        LOCAL_VALUES = Normalize_Lavalette(r,param)
        plt.plot(r,LOCAL_VALUES[0],linewidth=2)
        plt.scatter(r,LOCAL_VALUES[0],marker=shapes[PARÁMETROS.index(param)],label=str(round(param,2)))
    
    plt.legend(title=r'$b$', loc='best', bbox_to_anchor=(1, 0.525))
    plt.xlabel('Rango', fontproperties=label_font)
    plt.ylabel('Probabilidad', fontproperties=label_font)
    plt.title('Lavalette', fontproperties=label_font)
    

    
####### GERMIBETA #######
#def Germibeta(r,b,a):
#    N = max(r)
#    return [((N + 1 - i) / (i ** a)) ** (b) for i in r]

def Germibeta(r,b,a):
    N = max(r)
    return [((N + 1 - i) ** (b)) / (i ** a) for i in r]


def Normalize_Germibeta(r,b,a):
    Anormal_list = Germibeta(r,b,a)
    #print(Anormal_list)
    SUM = sum(Anormal_list)
    K = SUM ** (-1)

    
    return [[K * i for i in Anormal_list],b,a,K]


def Germi_plot_compare(b_1,a_1,b_2,a_2):
    #global r
    PARÁMETROS =  [[b_1,a_1],[b_2,a_2]]
    
    for param in PARÁMETROS:
        LOCAL_VALUES = Normalize_Germibeta(r,param[0],param[1])
        plt.plot(r,LOCAL_VALUES[0],linewidth=2)
        plt.scatter(r,LOCAL_VALUES[0],marker=shapes[PARÁMETROS.index(param)],
                    label=r'$b$ = {}  $a$ = {}'.format(str(round(param[0],2)),str(round(param[1],2))))
    
    plt.legend(title="Parámetros", loc='best', bbox_to_anchor=(1, 0.525))
    plt.xlabel('Rango', fontproperties=label_font)
    plt.ylabel('Probabilidad', fontproperties=label_font)
    plt.title('Germibeta', fontproperties=label_font)
    
    
    
    
    
####################### Para el ajuste de curva #############################
def Germi_Logaritmica_lmfit(x,K,b,N,a):
    return np.log(K) + (b * np.log(N+1-x)) - (a * np.log(x))


def Germi_Beta_lmfit(x,K,b,a,N):
    return K * (((N + 1 - x) ** b) / (x ** a))
    #return K * ((N + 1 - x) ** b)
    
    
r_prima = [i for i in range(1,513)]



def counter_natural(N,M):
    assert M > N, "El segundo número debe ser mayor al primero XD"
    K = N-1
    return int((M*(M+1)/2) - (K*(K+1)/2))

