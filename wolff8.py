import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import matplotlib.patches as mpatches
import random
from numpy.random import rand
import os
from PIL import Image
from datetime import datetime
import time
from scipy.optimize import curve_fit
from functools import cache
from tqdm import tqdm
from collections import Counter,deque


def distancia_corta(i,lista):
    """
    Calcula la distancia más corta entre un nodo i al resto de nodos, segun las conexiones dadas en lista
    
    Argumentos:
    i: int
        Nodo a partir del cual se quieren estudiar las distancias
    lista: array
        Lista de conexiones

    Output:
    distancias: dict
        Diccionario con keys un nodo j y como values la distancia de i a j

    """
    visitados=set([i])
    cola=deque([(i,0)]) #nodo actual // distancia
    distancias={}
    while cola:
        actual, distancia= cola.popleft()
        distancias[actual]=distancia
        for vecino in lista[actual]:
            if vecino not in visitados:
                visitados.add(vecino)
                cola.append((vecino,distancia+1))
                
    return distancias

def distancia_media(nodos,lista):
    """
    Cálculo de la distancia promedio entre nodos.
    Para cada nodo, se calcula la distancia al resto de ellos y se realiza la media.

    Argumentos:
    nodos: array
        Vector de nodos, donde en cada componente k se recoge el espín del nodo k
    lista: array
        Lista de conexiones
    
    Output:
    d: float
        Distancia promedio entre nodos.
    """
    lenv=len(nodos)
    suma=0
    ncon=0
    print("Calculo de la distancia promedio. ")
    for i in tqdm(range(lenv)):
        distancias=distancia_corta(i,lista)
        for j in range(i+1,lenv):
            suma+=distancias[j]
            ncon+=1
    d=suma/ncon
    return d


def plot_red(nodos, lista, L,directorio):
    """
    Visualización de la red. Crea una imagen para visualizar la red y sus conexiones
    
    Argumentos:
    nodos: array
        Vector de nodos, donde en cada componente k se recoge el espín del nodo k
    lista: array
        Lista de conexiones.
    L: int
        Longitud de la red.
    directorio: str
        Directorio donde guardar la imagen.
    """
    # Calculamos las posiciones (x,y) de cada nodo en la cuadrícula
    pos = {}
    for k in range(len(nodos)):
        fila = k // L
        col = k % L
        pos[k] = (col, -fila)  # -fila para que la gráfica tenga fila 0 arriba

    plt.figure(figsize=(8,8))
    
    # Dibujar las conexiones (aristas)
    for k, vecinos in enumerate(lista):
        x0, y0 = pos[k]
        for v in vecinos:
            x1, y1 = pos[v]
            # Dibujamos línea si k < v para evitar doble dibujo
            if k < v:
                plt.plot([x0, x1], [y0, y1], 'k-', lw=0.7, alpha=0.5)

    # Dibujar los nodos
    xs = [pos[k][0] for k in range(len(nodos))]
    ys = [pos[k][1] for k in range(len(nodos))]
    plt.scatter(xs, ys, c='red', s=40)

    plt.title("Visualización de la red LxL")
    plt.axis('off')  
    plt.gca().set_aspect('equal', adjustable='box')
    plt.savefig(directorio+"red.png",dpi=500)



def estandar(lenv,L,n,K=2):
    """
    Cálculo de las conexiones en el caso de primeros K/2 vecinos.

    Se basa en la transformación del nodo en forma vectorial al nodo en forma "matricial".

    Dado un nodo k, se construyen las componentes (i1,i2,...,in), donde n es la dimensionalidad de la red.
    Estas componentes se escogen de manera que i1*L^(n-1) + i2*L^(n-2) + ... + i_{n-1} * L + i_n = k
    Las conexiones establecidas serán con los nodos (i1 +- j, i2, ..., in), (i1,i2 +- j,...in), ... ,(i1, i2, ..., in +- j), donde
    j va desde 1 hasta K/2, teniendo en cuenta las condiciones periódicas.
    Una vez establecidas las conexiones, se transforman a forma vectorial sin más que calcular el polinomio anterior cambiando los 
    coeficientes correspondientes.

    Argumentos:
    lenv: int
        Longitud del vector de nodos.
    L: int
        Longitud de la red.
    n: int
        Dimensionalidad de la red
    K: int (par)
        Define el número de primeros vecinos de las conexiones: K/2.

    Output:
    lista:
        Lista de conexiones.
    
    """
    lista=[]
    for k in range(lenv):
        vcon=[]
        coeffs=[0]*n
        kcalc=k
        for i in range(n):
            coeffs[i]=kcalc//(L**(n-i-1))
            kcalc=kcalc%(L**(n-i-1))


        def calc(coeffs,ind,suma,n,L):
            conex=0
            coeffs[ind]=(coeffs[ind]+suma)%L
            for i in range(n):
                conex=conex+coeffs[i]*L**(n-i-1)
            return conex
        
        for i in range(n):
            for j in range(int(K/2)):
                vcon.append(calc(coeffs.copy(),i,j+1,n,L))
                vcon.append(calc(coeffs.copy(),i,-j-1,n,L))
        lista.append(vcon)
    return lista

def small_world(lenv,params):
    # Parámetros
    #N número de nodos
    #k cada nodo está conectado a k vecinos
    #p probabilidad de reconexión

    # Crear red small-world con el modelo de Watts-Strogatz
    N=lenv
    k,p=params
    lista=[]
    for i in range(N):
        conex=[]
        for j in range(1,int(k/2)+1):
            conex.append((i+j)%N)
            conex.append((i-j)%N)
        lista.append(conex)

    for i in range(N):
        conex=lista[i]
        for k in range(len(conex)):
            j=conex[k]
            if j>i:
                if np.random.rand() < p:
                    opciones=[i for i in range(N) if i not in conex]
                    conex[k]=random.choice(opciones)
        lista[i]=conex
    return lista


def small_world2(lenv,L,params):
    # Parámetros
    #N número de nodos
    #k cada nodo está conectado a k vecinos
    #p probabilidad de reconexión

    # Crear red small-world con el modelo de Watts-Strogatz
    """
    
    """
    N=lenv
    K,p=params
    lista=estandar(lenv,L,2,K)
    for i in tqdm(range(N)):
        conex=lista[i]
        for j in conex:
            if j>i:
                if np.random.rand() < p:
                    opciones=[l for l in range(N) if l not in conex]
                    elegido=random.choice(opciones)
                    lista[i].append(elegido)
                    lista[elegido].append(i)
                    lista[i].remove(j)
                    lista[j].remove(i)
    return lista
def small_world1d(lenv,L,params):
    # Parámetros
    #N número de nodos
    #k cada nodo está conectado a k vecinos
    #p probabilidad de reconexión

    # Crear red small-world con el modelo de Watts-Strogatz
    """
    
    """
    N=lenv
    K,p=params
    lista=estandar(lenv,L,1,K)
    for i in tqdm(range(N)):
        conex=lista[i]
        for j in conex:
            if j>i:
                if np.random.rand() < p:
                    opciones=[l for l in range(N) if l not in conex]
                    elegido=random.choice(opciones)
                    lista[i].append(elegido)
                    lista[elegido].append(i)
                    lista[i].remove(j)
                    lista[j].remove(i)
    return lista

def scale_free(lenv,L,params):
    lista=estandar(lenv,L,2,2)
    #lista=[[0]]*lenv
    m=params
    for i in range(lenv):
        vcon=[]
        for _ in range(m):
            grados=[len(conexiones) for conexiones in lista]
            grados = np.array(grados, dtype=float)
            probs = grados / grados.sum()
            vcon= np.random.choice(len(grados), size=m, replace=False, p=probs)
            vcon=vcon.tolist()
        lista[i].extend(vcon)
    grados=[len(conexiones) for conexiones in lista]
    print("MAXIMO DE CONEXIONES")
    print(max(grados))
    print(min(grados))


    # Aplanar la lista de listas
    aplanado = [nodo for sublista in lista for nodo in sublista]

    # Contar ocurrencias
    conteo = Counter(aplanado)

    # Ordenar por nodo (clave)
    nodos = sorted(conteo.keys())
    apariciones = [conteo[n] for n in nodos]

    # Plotear
    plt.figure(figsize=(10, 5))
    plt.bar(nodos, apariciones)
    plt.xlabel("Nodo")
    plt.ylabel("Número de apariciones")
    plt.title("Frecuencia de aparición de cada nodo")
    plt.grid(True, linestyle='--', alpha=0.3)
    plt.show()

    return lista



def conexiones(lenv,L,config):
    """
    Cálculo de las conexiones según el tipo de red.

    Además, devuelve un vector donde cada componente k contiene lista de las conexiones de cada nodo k.
    Estas conexiones vienen dadas dependiendo del tipo de red especificada.

    Argumentos:
    lenv: int
        Cantidad de nodos
    L: int
        Longitud de la red.
    config: array
        Vector de la forma [tipo, params], donde 
        tipo: str
            Determina el tipo de red: 
                Conexiones estándar: "std"
                Red Small World: "sw2"
                Red Free Scale: "sf"
        params: array
            Determina los parámetros según el tipo de red.
                "std": n
                    n: int
                        Dimensionalidad de la red
                "sw": [N,k,p]
                    N: int
                        Número de nodos
                    k: int (par)
                        Número de conexiones
                    p: float
                        Probabilidad de reconexión
                "sf": m
                    m: int
                        Número de conexiones
    
    Output:
    lista: vector
        Vector con las conexiones de cada nodo
    """
    tipo,params=config
  
    if tipo=="std":
        lista=estandar(lenv,L,params)
    elif tipo=="sw":
        lista=small_world(lenv,params)
    elif tipo=="sw2":
        lista=small_world2(lenv,L,params)
    elif tipo=="sw1d":
        lista=small_world1d(lenv,L,params)
    elif tipo=="sf":
        lista=scale_free(lenv,L,params)
    return lista 


def energia_vect(nodos,lista,J=1):
    """
    Calcula la energia en funcion del vector nodos.

    Argumentos:
    nodos: vector
        Vector de nodos, donde en cada componente k se recoge el espín del nodo k
    lista: array
        Lista de conexiones
    J: float
        Constante de acoplamiento en el modelo de Ising

    Devuelve:
    energia: float
        Energía calculada según el modelo de Ising
    """
    energia=0
    for i in range(len(lista)):
        for j in lista[i]:
            if j>i:
                energia+=-J*nodos[i]*nodos[j]
    return energia




def simulacion(nodos,lista,t,iters,directorio,J=1):
    """
    Ejecución del Algoritmo de Wolff

    Argumentos:
    nodos: array
        Vector de nodos, donde en cada componente k se recoge el espín del nodo k
    lista: array
        Lista de conexiones.
    t: float
        Temperatura    
    iters: int
        Número de iteraciones del algoritmo de Wolff
    directorio: str
        Directorio donde guardar los datos de la ejecución.
    J: float
        Constante de acoplamiento en el modelo de Ising
    
    Devuelve:
    magnetizacion: vector
        Magnetizacion promedio del sistema en cada iteracion
    energias: vector
        Energía promedio del sistema en cada iteracion

    """
    print("Comienzo de simulacion")

    beta=1/t
    magnetizacion=np.zeros(iters)
    energias=np.zeros(iters)

    E=energia_vect(nodos,lista)
    
    #ALGORITMO DE WOLFF 
    for m in tqdm(range(iters)):
        k0=random.randint(0,len(nodos)-1) #Selección de nodo aleatorio
        estado=nodos[k0]
        cluster=[k0]
        visited=np.zeros_like(nodos,dtype=bool)
        stack=[k0]
        visited[k0]=True
        prob=1-np.exp(-2*beta*J)

        while stack:
            # con esta funcion se estudian los nodos aun no estudiados.
            # Se estudian las conexiones de los nodos añadidos en la iteración anterior.
            # En caso de no añadirse ningun nodo nuevo, se detiene.
            k=stack.pop()
            nodos[k]
            for conex in lista[k]: 
                if not visited[conex] and nodos[conex]==estado:
                    if np.random.rand()<prob:
                        cluster.append(conex)
                        stack.append(conex)
                        visited[conex]=True
           
        indices=np.array(cluster)
        
        deltaE=0
        for k in indices:
            for conex in lista[k]:
                if not visited[conex]:
                    deltaE += 2*J*estado*nodos[conex]
        E=E+deltaE
        nodos[indices]=nodos[indices]*-1 # se dan la vuelta a los nodos del cluster
        energias[m]=E
        magnetizacion[m]=np.abs(np.mean(nodos))
    
    if os.path.exists(directorio+"temps.npy"):
        temps=np.load(directorio+"temps.npy")
        temps=np.append(temps,t)
        np.save(directorio+"temps.npy",temps)
    else:
        temps=np.array([t])
        np.save(directorio+"temps.npy",temps)
    ind=len(temps)-1
    np.save(directorio+"energias"+str(ind)+".npy",energias)
    np.save(directorio+"mags"+str(ind)+".npy",magnetizacion)

    return magnetizacion,energias#,mat_nodos



def programa(L,iters,temps,directorio_fecha,directorio_ejec,config,J=1):
    """
    Ejecuta la simulacion del modelo de Ising utilizando el algoritmo de Wolff.
    En una carpeta guarda los datos y gráficas de interes del sistema.

    Argumentos:
    filas: int
        Número de filas de la red
    cols: int
        Número de columnas de la red
    iters: int
        Número de iteraciones del algoritmo de Wolff
    temps: vector
        Vector con las temperaturas a estudiar
    directorio_fecha: str
        Directorio donde se almacenan los datos y gráficas de la simulacion para unos parametros determinados
    directorio_ejec:
        Subcarpeta de directorio_fecha donde se guardan los datos y gráficas de cada ejecucion

    Devuelve:
    mat_e: matriz
        Matriz donde en cada fila i se guardan las energias en cada iteración para la temperatura temp[i]
    mat_m: matriz
        Matriz donde en cada fila i se guardan las magnetizaciones en cada iteracion para la temperatura temp[i]
    matrices: lista
        Lista donde en cada componente i se guarda una matriz. En esta matriz para la temperatura temp[i],
        en cada fila k se guarda el vector nodos en la iteracion k
    """
    
    directorio=directorio_fecha+directorio_ejec+"/"
    os.makedirs(directorio,exist_ok=False) #carpeta donde guardar las imagenes
    plt.close()

    tipo,params=config
    if tipo=="std":
        n=params
        #Creación del vector de nodos
    if tipo=="sw1d":
        n=1
    else:
        n=2
    nodos=np.ones(L**n)
    lenv=L**n

    #Creación de las conexiones de la red
    if os.path.isfile(directorio_fecha+"lista.npy"):
        lista=np.load(directorio_fecha+"lista.npy", allow_pickle=True)
        print("Lista ya existente")
    else:
        print("Lista creada por primera vez")

        lista=conexiones(lenv,L,config) #lista de conexiones`
        np.save(directorio_fecha+"lista.npy",np.array(lista, dtype=object))
        plot_red(nodos,lista,L,directorio_fecha)

    #d=distancia_media(nodos,lista)
    #print(d)
    #Archivo de texto donde se guardan datos generales
    f = open(directorio+"datos.txt","w")
    f.write("Numero de iteraciones: "+str(iters)+"\n")
    f.write("J="+str(J)+"\n")
    f.write("Dimension de la red: "+str(L)+"x"+str(L)+"\n")
    f.write("Filas:"+str(L)+"\n")
    f.write("Temperaturas: "+str(temps)+"\n")
    f.write("Tipo de red. n="+str(n))
    f.close()
    
    #Simulación
    mat_e=np.zeros((len(temps),iters)) #en cada fila, el vector energias para la temperatura t
    mat_m=np.zeros((len(temps),iters)) #en cada fila, el vector magnetizaciones para la temperatura t
    for k in range(len(temps)):
        t=temps[k]
        print(str(k)+"/"+str(len(temps))+". Temperatura: ", str(round(t,4)))
        mat_m[k,:],mat_e[k,:]=simulacion(nodos.copy(),lista,t,iters,directorio,J) 
           
    
    print("Creacion de plots")
    tipo,params=config
    if tipo=="std":
        n=params
        titulo=str(n)+"D2"
    elif tipo=="sw":
        k,p=params
        titulo="k="+str(k)+", p="+str(p)
    elif tipo=="sw2":
        k,p=params
        titulo="k="+str(k)+", p="+str(p)
    elif tipo=="sf":
        m=params
        titulo="m="+str(m)



    """ for i in tqdm(range(len(temps))):
        t=temps[i]
        plt.plot(mat_e[i,:])
        plt.title("E vs iters, T="+str(round(t,3))+" L= "+str(L)+". "+titulo)
        plt.xlabel("iteracion")
        plt.ylabel("energia")
        plt.savefig(directorio +str(i)+"energias.png",dpi=50)
        plt.close()
        plt.plot(mat_m[i,:])
        plt.title("M vs iters T="+str(round(t,3))+" L= "+str(L)+titulo)
        plt.xlabel("iteracion")
        plt.ylabel("Magnetizacion")
        plt.savefig(directorio +str(i)+"mags.png",dpi=50)
        plt.close() """
    return mat_e, mat_m

    



os.chdir(r"C:\Users\panti\OneDrive\Documentos\Pablo\Universidad\2024-25\Semestre2\TFG Física")
######################################
#########---CONFIGURACION---##########
######################################

#Dimensión de la red
L=200
print("L="+str(L))
#Constante J
J=1

####Tipo de red
##ESTANDAR
n=2
params_std=n

##Small World
k=4
p=0
params_sw=[k,p]

##Scale Free
m=20
params_sf=m

#Elección de tipo
tipo="sw2"


#Selección de temperaturas e iteraciones

lim_inf=3#float(input("Limite inferior? ")) #2.2  #2.25 #1.5 #2.8  #0.01   ##3D 3
lim_sup=5#float(input("Limite superior? ")) #2.8  #2.35 #2.2 #3.5  #3           5
cant=101#int(input("Cantidad? "))           #50   #20  #25  #25   #50           101
iters=20000#int(input("Iteraciones? "))

######################
######################
######################

#Creación de carpeta
if tipo=="std":
    params=params_std
    nombre_carpeta=str(n)+"d"
    carpeta_config="ZZZ_"+nombre_carpeta+"/wolff_"+str(L)+"f"

elif tipo=="sw":
    params=params_sw
    nombre_carpeta=str(k)+"k_"+str(p).replace(".","_")+"p"
    carpeta_config="ZZZ_"+nombre_carpeta+"/wolff_"+str(L)+"f"

elif tipo=="sw2":
    params=params_sw
    nombre_carpeta="SW2"#str(k)+"k_"+str(p).replace(".","_")+"pSW2"
    carpeta_config="ZZZ_"+nombre_carpeta+"/wolff_"+str(k)+"k_"+str(p).replace(".","_")+"p"
elif tipo=="sw1d":
    params=params_sw
    nombre_carpeta="SW1d"#str(k)+"k_"+str(p).replace(".","_")+"pSW2"
    carpeta_config="ZZZ_"+nombre_carpeta+"/wolff_"+str(k)+"k_"+str(p).replace(".","_")+"p"
elif tipo=="sf":
    params=params_sf
    nombre_carpeta=str(m)+"m_SF"
    carpeta_config="ZZZ_"+nombre_carpeta+"/wolff_"+str(L)+"f"




config=[tipo,params]

directorio_config="fotos/"+carpeta_config+"/"
os.makedirs(directorio_config,exist_ok=True) 



temps=np.linspace(lim_inf,lim_sup,cant)
temps=temps[(temps>4.5) & (temps<7)]

if n==3:
    temps=np.arange(3,6.02,0.02)
    temps=np.append(temps,[4.505,4.51,4.515])
    temps=np.sort(temps)
    lim0=4.3
    lim1=4.49
    iters=1200
    temps=temps[(temps>=lim0) & (temps<=lim1)]
if n==4:
    temps=np.arange(6,7.52,0.02)
    temps=np.append(temps,[6.665,6.675,6.685])
    temps=np.sort(temps)
    lim0=5.8
    lim1=6.42
    iters=1300
    temps=temps[(temps>=lim0) & (temps<=lim1)]
if tipo=="sw2":
    #if p==0.4 and k==2:
    temps=np.arange(2,8.02,0.02)
    #temps=np.append(temps,[])
    temps=np.sort(temps)
    lim0=6.1
    lim1=8
    iters=50000
    temps=temps[(temps>=lim0) & (temps<=lim1)]

if tipo=="sw1d":
    temps=np.arange(0,8.02,0.02)
    #temps=np.append(temps,[])
    temps=np.sort(temps)
    lim0=6.1
    lim1=8
    iters=70000
    temps=temps[(temps>=lim0) & (temps<=lim1)]


#temps=temps[len(temps)-2:]
np.save(directorio_config+"filas.npy",L)
np.save(directorio_config+"tipo.npy",tipo)
np.save(directorio_config+"params.npy",params)

ruta = directorio_config
carpetas = [nombre for nombre in os.listdir(ruta)
            if os.path.isdir(os.path.join(ruta, nombre))]
carpetas = [int(numeric_string) for numeric_string in carpetas if numeric_string.isdigit()]
if len(carpetas)==0:
    directorio_ejec=str(1)
else:
    directorio_ejec=str(max(carpetas)+1)
print(directorio_ejec)
print("n="+str(n))
mat_e,mat_m=programa(L,iters,temps,directorio_config,directorio_ejec,config,J)
