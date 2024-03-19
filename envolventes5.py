# Graficas de las cotas, envolventes y todos los puntos de los datos experimentales

import numpy as np
import matplotlib.pyplot as plt

# Directorio de donde se cargan datos y se va a guardar la imagen
# directorio_destino = r"C:\Users\eduar\Desktop\Ejercicios_Python\porosidades_max_y_min\HJ1468_M_total"
directorio_destino = r'C:\Users\eduar\Desktop\Ejercicios_Python\porosidades_2D\HJ1414_M'

##################################################################
#                     Datos experimentales
##################################################################
poro_exp = 0.23 # Porosidad experimental de la muestra cilindrica

# Porosidad tendencia de los datos (porosidad del prisma)
m = 0.246182

resolucion = 15

limite_grafica = 1.0
##################################################################

datos = np.loadtxt(directorio_destino + "\porosidades2D.txt")

serie_poro = datos/100.0

num_experimentos = len(serie_poro[1])

###############################################################################

# Obtener los valores mínimo y máximo de cada columna (eje 0)
max_values = np.max(serie_poro, axis=0)
min_values = np.min(serie_poro, axis=0)

num_submuestras = len(max_values)

##################################################################

# Unidades

micra_cub = 1e-12 # cm3
# mm_cub = 1e13 # micras cubicas
cm_cub = 1e12 # micras cubicas
micra = 1e-4 # cm
milimetro = 1000 # micras
cm = 10000 # micras

# Calculamos los valores del eje x
long = resolucion*micra*2*num_submuestras
paso = resolucion*micra*2

vector_x = np.arange(paso, long+paso, paso)

m_max = m

m_min = m

# Datos experimentales para la cota superior 

# Buscamos el índice donde las porosidades son diferente de cero
index_max = next((i for i, x in enumerate(max_values) if x <= 0.99), None)

# Calculamos el tamaño de poro máximo (tau = X0)
tau_max = index_max*paso # cm

# Coordenadas del punto de referencia
Lmax_end = vector_x[-1]
Pmax_end = max_values[-1]

# Buscamos el índice donde las porosidades son diferente de cero
index_min = next((i for i, x in enumerate(min_values) if x != 0), None)

# Calculamos el tamaño de poro mínimo (tau = X0)
tau_min = index_min*paso # cm

# Coordenadas del punto de referencia
Lmin_end = 1.0#vector_x[-1]
Pmin_end = 0.21#min_values[-1] 

##################################################################
# Definimos las funciones de las envolventes
##################################################################

def upper_bound(L, k, m, tau):
    if L <= tau:
        return 1
    else:
        return m + ((1.0-m)*np.exp(-k*(L-tau)))
    
def lower_bound(L, k, m, tau):
    if L <= tau:
        return 0
    else:
        return  m - (m*np.exp(-k*(L-tau)))

# Definimos las funciones de las constantes de proporcionalidad de las cotas
def kmax(L, P, m, tau):
    return -np.log(abs((P-m)/(1.0-m))) / (L-tau)

def kmin(L, P, m, tau):
    return -np.log(abs((P-m)/-m)) / (L-tau)
    
###############################################################################
# Constante de proporcionalidad de las cota superior
k_max = kmax(Lmax_end, Pmax_end, m_max, tau_max)

# Constante de proporcionalidad de la cota inferior
k_min = kmin(Lmin_end, Pmin_end, m_min, tau_min)

##################################################################
#                       G r a f i c a s
##################################################################
# Crear la figura y los ejes
fig, ax = plt.subplots(figsize=(4.0, 3.0))

##################################################################
# Puntos de los datos experimentales
##################################################################
# Crear un rango de colores utilizando el colormap 'jet'
num_lines = len(serie_poro)
colores = plt.cm.jet(np.linspace(0, 1, num_lines))

# Graficar los datos con colores del colormap
for i in range(num_lines):
    ax.plot(vector_x, serie_poro[i], '-', color=colores[i], linewidth = 0.1)

###############################################################################    
# Bounds
#############################################################################
# Crear un rango de valores de x para graficar las cotas superior e inferior
t_values = np.arange(paso, limite_grafica, paso)

#############################################################################
# Upper bound
#############################################################################
# Calcular los valores de la cota superior para cada valor de L
cota_superior = [upper_bound(L, k_max, m_max, tau_max) for L in t_values]
# Graficar la cota superior
ax.plot(t_values, cota_superior, '-', color='black' , linewidth = 1.25, label = 'Upper Bound')

#############################################################################
# Lower bound
#############################################################################
# Calcular los valores de la cota inferior para cada valor de L
cota_inferior = [lower_bound(L, k_min, m_min, tau_min) for L in t_values]
# Graficar la cota inferior
plt.plot(t_values, cota_inferior, '-', color = 'grey' , linewidth = 1.25, label = 'Lower Bound')
#############################################################################
# Crear el gráfico de la porosidad computacional
grafica_poro_comp = [m for L in t_values]
plt.plot(t_values, grafica_poro_comp, '--', color = 'magenta' , linewidth = 1.0, label = '$m = \phi_{\mu CT}$')

#############################################################################
# Crear el gráfico de la porosidad experimental
grafica_poro_exp = [poro_exp for L in t_values]
plt.plot(t_values, grafica_poro_exp, '-.', color = 'green' , linewidth = 1.0, label = '$\phi_{exp}$')

###############################################################################
# Agregar un colorbar a la derecha
sm = plt.cm.ScalarMappable(cmap='jet', norm = plt.Normalize(vmin=1, vmax=num_lines))
sm.set_array([])
cbaxes = fig.add_axes([0.91, 0.1, 0.03, 0.8])
cbar3 = plt.colorbar(sm, cax = cbaxes, location = 'right', pad=0.01)
cbar3.set_label('Sequence of experiments', fontsize=10)
cbar3.ax.tick_params(labelsize=8)  # Tamaño de la fuente del colorbar
cbar3.ax.invert_yaxis()  # Si quieres que los ticks del colorbar vayan de arriba a abajo
###############################################################################
# Configurar etiquetas y título
ax.set_xlabel('Scale [cm]', fontsize=10)
ax.set_ylabel('Porosity', fontsize=10)
ax.set_title('Sample 2 (Resolution: 15 $\mu$m/voxel)', fontsize=10)
ax.legend(fontsize=8)
# ax.grid(True)
# Ajustar el tamaño de las etiquetas y números de los ejes
ax.tick_params(axis='both', labelsize=8)
ax.tick_params(axis='both', which='major', labelsize=10)

# Guardar la grafica en una carpeta 
plt.savefig(directorio_destino + '\grafica_envolventes.png', bbox_inches='tight', dpi=300)

print('kmax = ', round(k_max, 2))
print('kmin = ', round(k_min, 2))
print('m = ', round(m, 2))
print('Tau_max = ', round(tau_max, 4))
print('Tau_min = ', round(tau_min, 4))
print('(Lmax,Pmax)=', '(', Lmax_end, ',', round(Pmax_end, 2), ')')
print('(Lmin,Pmin)=', '(', Lmin_end, ',', round(Pmin_end, 2), ')')