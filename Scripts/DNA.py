# Medición de tiempo de procesamiento y request para extracción de datos desde URL
from tqdm import tqdm
import requests, os
from lxml import html
import pandas as pd

# Descompresión de archivos
import gzip
import shutil

# Bio: Bioinformatics with python
import Bio
#from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import SimpleFastaParser
from Bio import SeqIO

from collections import Counter
import json



# Variables para la descarga
DATA_FOLDER = "../Data"
DATA_SUBFOLDER = os.path.join(DATA_FOLDER,"Chromosomes")
EXISTS_ALREADY = False

### Variables para la traducción
DICS = {'Enlace': {'A': 0, 'T': 0, 'C': 1, 'G': 1},
        'P&P' : {'A': 0, 'T': 1, 'C': 1, 'G': 0}}

### Listas con valores para limpieza
LETTER_ROUGH = ['A','a','T','t','C','c','G','g','N','OTHER']
LETTERS_CLEANED = ['A','T','C','G','N','OTHER']

def download(url, fname):
    """
    Función para descargar la información desde una liga. Como input requiere:
        url: de dónde descargar
        fname: Filename. Nombre del archivo dentro del sistema. Por defecto se descarga
            en la misma carpeta donde esté el notebook.
    """
    resp = requests.get(url, stream=True)
    total = int(resp.headers.get('content-length', 0))
    with open(fname, 'wb') as file, tqdm(
            desc=fname,
            total=total,
            unit='iB',
            unit_scale=True,
            unit_divisor=1024,
    ) as bar:
        for data in resp.iter_content(chunk_size=1024):
            size = file.write(data)
            bar.update(size)
            
            
def download_dna(DATA_FOLDER, DATA_SUBFOLDER):
    """
    Descarga la versión más reciente del genoma humano desde la página del NCBI
    y lo guarda en la carpeta de 'Data' del proyecto. No se incluyen estos datos
    en el repositorio porque el añadirlos lo haría innecesariamente pesado.
    """
    global EXISTS_ALREADY
    # Creamos una carpeta para guardar la información. En caso de existir
    # revisamos que no tenga ya un archivo con las características del genoma
    try:
        os.mkdir(os.path.join(DATA_FOLDER))
        os.mkdir(os.path.join(DATA_FOLDER,DATA_SUBFOLDER))
    except:
        for file in os.listdir(os.path.join(DATA_FOLDER)):
            if "GCA" in file:
                EXISTS_ALREADY = True
    
    if not EXISTS_ALREADY:
        # Nos dijirimos a la página de ncbi con la versión del cromosoma humano más reciente
        # no es necesario entrar en el detalle de estas funciones aunque se incluyen breves descripciones

        # Página base de dónde extraemos la información
        FTP_HOMOSAPIENS_ROOT = "https://ftp.ncbi.nih.gov/genomes/genbank/vertebrate_mammalian/Homo_sapiens/reference/"

        # Se visita la página con un request.get para indagar sobre el contenido
        # Esto se hace porque la información está en cambio constante y cada versión
        # tiene una nomenclatura diferente a las anteriores.
        HS_ROOT = requests.get(FTP_HOMOSAPIENS_ROOT)
        webpage = html.fromstring(HS_ROOT.content)

        # Tomamos la última versión y visitamos esa página
        LIST_CONTENTS = webpage.xpath('//a/@href')
        HYPER_LINK = LIST_CONTENTS[1]

        HS_ROOT = requests.get(FTP_HOMOSAPIENS_ROOT + HYPER_LINK)
        webpage_1 = html.fromstring(HS_ROOT.content)

        LIST_CONTENTS_1 = webpage_1.xpath('//a/@href')

        # Nos dirigimos al archivo con el formato fasta y comprimido como gz
        FASTA_FILE_NAME = HYPER_LINK[0:-1] + '_genomic.fna.gz'
        HUMAN_GENOME_FASTA = FTP_HOMOSAPIENS_ROOT + HYPER_LINK + FASTA_FILE_NAME

        # Descargamos la información mostrando una barra de progreso para que el usuario tenga conocimiento
        # de lo que está ocurriendo
        print("Descargando la versión más reciente del genoma humano, en formato fasta desde {}".format(FTP_HOMOSAPIENS_ROOT + HYPER_LINK))
        download(HUMAN_GENOME_FASTA, os.path.join(DATA_FOLDER,FASTA_FILE_NAME))
        
    else:
        # Durante las pruebas, se corrió esta sección del código en varias ocasiones y no siempre queremos
        # que se descarge.
        # También podría ocurrir que el usuario olvide cambiar el valor de la variable de descarga así que
        # prevenimos descargas innecesarias.
        print("Ya existe un archivo en la carpeta de destino. Por favor asegúrate que no tengas ya los datos")    
        
        
    ######################################################################
    #### Ahora extraemos los datos y los guardamos como archivos .txt ####
    ######################################################################
    
def unzip_dna(DATA_FOLDER):
    FASTA_FILE_NAME = [file for file in os.listdir(DATA_FOLDER) if ".gz" in file][0]
    with gzip.open(os.path.join(DATA_FOLDER,FASTA_FILE_NAME), "rt") as handle:
        print(sum(len(r) for r in SeqIO.parse(handle, "fasta")))
        
        
def unzip_dna_store(DATA_FOLDER):
    """
    Descompresión del genoma humano que previamente guardamos en la carpeta 'Data'
    No es necesario que se haya descargado con las funciones incluídas en este proyecto
    pero sí que tenga una extensión .gz y que esté en la carpeta correcta.
    
    Para evitar una descompresión innecesaria se verifica primero que no haya archivos descomprimidos.
    """
    FASTA_FILE_NAME = [file for file in os.listdir(DATA_FOLDER) if ".gz" in file][0]
    TARGET_NAME = FASTA_FILE_NAME[0:-3] + ".txt"
    
    # Verificamos que no haya un archivo descomprimido en la carpeta
    if (TARGET_NAME in os.listdir(DATA_FOLDER)):
        print("La carpeta seleccionada ya tiene un archivo descomprimido con las características necesarias. Favor de verificar")
        print(TARGET_NAME)
    else:
        with gzip.open(os.path.join(DATA_FOLDER,FASTA_FILE_NAME), 'rb') as f_in:
            with open(os.path.join(DATA_FOLDER,TARGET_NAME), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
                
                

#######################
    
    
    
# Conversión de datos binarios a decimales
def binary_decimal(input_num):
    """
    Función que recibe un número expresado en base binaria y lo convierte a base decimal
    """
    assert type(input_num) == str, "El input debe ser un número binario dado como un string"
    b_num = list(input_num)
    value = 0
    
    for i in range(len(b_num)):
        digit = b_num.pop()
        if digit == '1':
            value = value + pow(2, i)
            
    return value

# Mapping of a string to another through the application of a dictionary to all the values (letters)
def map_string_dict(string,diccionary):
    """
    Mapeo de un string a otro mediante un diccionario.
    ABCD -> Por medio del diccionario Dic = {A:1, B:2, C:3, D:4} -> Regresa str(1234)
    """
    return "".join([str(diccionary.get(letter,"0")) for letter in string])


def sort_dict_values(DIC,by='values',rev=True):
    """
    Ordenamiento de los valores de un diccionario. 
    Por defecto se ordenan de acuerdo con los valores (parte derecha de la adjudicación)
    aunque esto puede cambiarse como argumento de la función.
    Únicamente habría que añadir el argumento 
        by='keys'
    """
    if by == 'values':
        return {k: v for k, v in sorted(DIC.items(), key=lambda item: item[1], reverse=rev)}
    elif by.lower() == 'keys':
        return {k: v for k, v in sorted(DIC.items(), key=lambda item: item[0], reverse=rev)}

    
# Definition of k-mers of variable length from a string
def string_kmers_parse(string,mapping_rule,k=3,deletion_rule='DNA'):
    # Since no ambiguity might arise due the use of other dictionaries we delete the undefined 
    # letters
    #string = string.replace(DELETE_ALPHABET[deletion_rule],"")
    
    #dicctionary = {}
    lista_valores = []
    
    # Quitamos el residuo para que se itere sin errores
    Residual = (len(string) % k)
    
    if Residual == 0:
        STRING = string
    else:
        STRING = string[0:-(len(string) % k)]
    
    #print(string,STRING)
    #print(STRING)
    #print(len(STRING))
    #print(len(STRING)/k)
    
    for iteracion in range(0,int(len(STRING)/k)):
        lista_valores.append(mapping_rule(STRING[0+(iteracion*k):k+(iteracion*k)]))
        #print(STRING[0+(iteracion*k):k+(iteracion*k)])
        #print(STRING)
        
    return lista_valores ,len(string) % k, sort_dict_values(dict(Counter(lista_valores)))


def random_kmer_string(String,mapping_rule,Lista):
    """
    Creamos ahora un arreglo de longitudes variables de un conjunto de posibilidades
    """
    lista_valores = []
    Filled = False
    Remainder = len(String)
    STRING = String[0:]
    
    while not Filled:
        Length = np.random.choice(Lista)
        #print(Length,STRING)
    
        if Length >= Remainder:
            lista_valores.append(mapping_rule(STRING))
            Filled = True
            
        else:
            lista_valores.append(mapping_rule(STRING[0:Length+1]))
            STRING = STRING[Length:]
            
        Remainder -= Length
            
    return lista_valores, sort_dict_values(dict(Counter(lista_valores)))
        
    
def fasta_genoma_HS_diccionario_cromosomas():
    
    """
    Recibe un archivo fasta previamente interpretado como un texto, aquí llamado 
    
    GCA_000001405.15_GRCh38_genomic.fna.txt
    
    Y genera un diccionario anidado con la estructura
    
    CHROMOSOME_DICTIONARY = {
        'cromosoma 1':{
        'cadena':[ ... ],
        'length': ##
        }
        ,
        'cromosoma 2':{
        'cadena':[ ... ],
        'length': ##
        }
    ...
    
        'cromosoma N':{
        'cadena':[ ... ],
        'length': ##
        }
    }
    """

    CHROMOSOME_DICTIONARY = {"chromosome 21":{'cadena':str(),'longitud':0}}
    
    with open("GCA_000001405.15_GRCh38_genomic.fna.txt") as handle:
        for title, seq in SimpleFastaParser(handle):

            if "chromosome 21" in title:
                # Cadena con el alias chromosome ## número de cromosoma
                SPLIT = str(title.split('chromosome ')[1]).replace(',','')
                
                try:
                    NUMBER = [int(i) for i in SPLIT.split() if i.isdigit()]
                    STRING = "chromosome " + str(NUMBER[0])
                except:
                    ST = [i for i in SPLIT.split()][0]
                    STRING = "chromosome " + ST

                # Crea un diccionario para guardar valores
                if not STRING in CHROMOSOME_DICTIONARY.keys():
                    CHROMOSOME_DICTIONARY[STRING] = {'cadena':str(),'longitud':0}

                # En caso de no requerir crear un diccionario; únicamente lo llenamos
                Secuencia = seq.upper().replace('N','')
                CHROMOSOME_DICTIONARY[STRING]['cadena'] += seq
                CHROMOSOME_DICTIONARY[STRING]['longitud'] += len(seq)            

    return CHROMOSOME_DICTIONARY


def store_dic_json(dictionary,name):
    """
    Guarda un diccionario como archivo .txt con una previa transformación a json
    esto para que se guarde como archivo 'serializable'
    
    Recibe dos argumentos:
        
        1.- El diccionario que quiere convertirse y guardarse.
        2.- El nombre bajo el cual se quiere guardar.
    """
    
    global Path_to_Sequence
    #name = variable_name(dictionary,globals())[0]
    with open(os.path.join(Path_to_Sequence, name + '.txt'), 'w') as file:
        file.write(json.dumps(dictionary))
        file.close()
        
def store_dic_json_path(path,dictionary,name):
    """
    Guarda un diccionario como archivo .txt con una previa transformación a json
    esto para que se guarde como archivo 'serializable'
    
    Recibe dos argumentos:
        
        1.- El diccionario que quiere convertirse y guardarse.
        2.- El nombre bajo el cual se quiere guardar.
    """
    
    #global Path_to_Sequence
    #name = variable_name(dictionary,globals())[0]
    with open(os.path.join(path, name + '.txt'), 'w') as file:
        file.write(json.dumps(dictionary))
        file.close()
        
        
def dic_whole_storing_json(DIC):
    for key in DIC.keys():
        store_dic_json(DIC[key],key)
        
        
def read_json_file(json_file):
    """
    Carga un archivo .txt con formato json
    Considerar por favor que nuevamente se asume que el archivo json_file.txt está guardado
    en la misma ruta donde se han estado guardando las secuencias
    """
    
    #name = variable_name(json_file,globals())[0]
    with open(os.path.join(Path_to_Sequence, json_file + '.txt')) as f:
        return json.load(f)
    
def read_json_file_path(Path_to_Sequence,json_file):
    """
    Carga un archivo .txt con formato json
    Considerar por favor que nuevamente se asume que el archivo json_file.txt está guardado
    en la misma ruta donde se han estado guardando las secuencias
    """
    
    #name = variable_name(json_file,globals())[0]
    with open(os.path.join(Path_to_Sequence, json_file + '.txt')) as f:
        return json.load(f)
    
def variable_name(obj, namespace):
    
    """
    Función que genera dos variables, la segunda es el nombre del objeto para después usarlo en 
    funciones que guarden los nombres
    """
    
    return [name for name in namespace if namespace[name] is obj]

def save_string_txt(PATH,string):
    """
    Stores a desired string to a .txt file.
    
    PLEASE CONSIDER: These are stored in D:\Documentos\Thesis\Bioinformatics_RobertoBastida\Dynamical Systems DNA\Data_Human_DNA\DNA_Extractions
    change this to any desired folder in the initial definition of variables in the current script
    """
    name = variable_name(string,globals())[0]
    
    #global Path_to_Sequence
    text_file = open(os.path.join(PATH, name + '.txt'), "wt")
    n = text_file.write(string)
    text_file.close()
    
    
def save_string_txt_name(PATH,name,string):
    """
    Stores a desired string to a .txt file.
    La diferencia con save_string_txt es que ésta función permite al usuario adjudicarle el nombre con el que desee
    guardar el texto
    """
    text_file = open(os.path.join(PATH, name + '.txt'), "wt")
    n = text_file.write(string)
    text_file.close()
    

def string_kmers_parse_tracker(string,mapping_rule,k=3):
    #dicctionary = {}
    lista_valores = []
    
    
    # Quitamos el residuo para que se itere sin errores
    STRING = string[0:-(len(string) % k)]
    #print(STRING)
    #print(len(STRING))
    #print(len(STRING)/k)
    
    for iteracion in range(0,int(len(STRING)/k)):
        lista_valores.append(mapping_rule(STRING[0+(iteracion*k):k+(iteracion*k)]))
        #print(STRING[0+(iteracion*k):k+(iteracion*k)])
        
    return lista_valores

def Display_First_N_Dictionary(DIC,N,columnas=['Keys','Values']):
    """
    Muestra las primeras N filas de un diccionario. Únicamente válido en notebooks 
    """
    Smaller = dict(zip(list(DIC.keys())[0:N],list(DIC.values())[0:N]))
    return pd.DataFrame(list(Smaller.items()),columns=columnas)

def list_to_string(lista):
    """
    Convierte una lista a string
    """
    s = str()
    
    for i in lista:
        s += str(i)
        #print(s)
        
    return s


def Convert(string):
    """
    Conversión de cadena de caracteres a lista separada por caracter
    """
    list1=[] 
    list1[:0]=string 
    return list1 


LISTA_LIMPIA = ['A','T','G','C']
DICT_LIMPIEZA = {'a':'A','t':'T','c':'C','g':'G'}
def Cleansing(String):
    """
    Como se menciona en repetidas ocasiones. Queremos remover las minúsculas así como caracteres
    no deseados (en ocasiones hay otras letras diferentes a las que consideramos en los        
    diccionarios).
    
    Si bien sería más sencillo (y rápido) aplicar el método .upper a la cadena de caracteres
    previamente depuradas de la letra N; eso no basta para identificar casos anómalos
    que potencialmente pueden alterar toda la ejecución
    """
    RESULTADO = []
    for letra in String:
        if letra in LISTA_LIMPIA:
            RESULTADO.append(letra)
        elif letra in list(DICT_LIMPIEZA.keys()):
            RESULTADO.append(DICT_LIMPIEZA[letra])
        else:
            pass
        
    return list_to_string(RESULTADO)
            

def Converter_and_Deletion(string):
    """
    Manera rápida de eliminar las letras N del conjunto
    """
    return(string.upper().replace("N",""))

def transformation_rule_(Rule,string):
    s = str()
    lista = list(string)
    counter_error = 0
    
    for j in lista:
        #print(j)
        #s += str(Rule[str(j.upper())])
        try:
            s += str(Rule[str(j.upper())])
        except:
            counter_error += 1
            pass
    
    print(counter_error)
    return s

    
def load_whole_chromosome():
    """
    Función que carga a memoria un archivo .txt
    """
    global Path_Sample_Sequence
    return open(os.path.join(Path_Sample_Sequence,"chromosome_21_.txt"),"r").read()

def load_txt_file(PATH,file):
    """
    Función que carga a memoria un archivo .txt
    """
    #global Path_Sample_Sequence
    return open(os.path.join(PATH,file),"r").read()