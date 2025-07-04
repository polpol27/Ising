# Ising
Estudio del modelo de Ising mediante el algoritmo de Wolff, realizado para el TFG de Física en la USC.

Para poder utilizar este programa, se especifican los diferentes parámetros al final del mismo:
L: longitud de la red
J: constante de acoplamiento

Red interacción primeros vecinos:
n: dimensión de la red
tipo: "std"

Red Small World 1D:
k: la red se conecta inicialmente a los k/2 primeros vecinos a cada lado
p: probabilidad de reconexión
tipo: "sw1d" Red 1D
tipo: "sw2" Red 2D

Temperaturas:
Se guardan en el vector "temps"
lim_inf: límite inferior
lim_sup: límite superior
cant: cantidad de temperaturas a analizar entre lim_inf y lim_sup

iters: número de iteraciones

El programa creará una carpeta en el directorio actual con un nombre específico, según la configuración realizada. En su interior, se guardarán los datos característicos del modo.
En cada ejecución, este creará una carpeta numérica donde se guardarán los datos de la magnetización y de la energía en cada iteración, para cada temperatura.
Opcionalmente, se guardan gráficas con estas magnitudes frente a las iteraciones, para observar la llegada al equilibrio.
