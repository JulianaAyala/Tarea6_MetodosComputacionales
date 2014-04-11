
###Autores:
###Juliana M. Ayala
###David Aleman
## codigo desarrollado en 2013-2
###PYTHON!
#########################################################################

import numpy as np
import sys, string, os
from StringIO import StringIO

#defino una funcion que cada vez que le meta el nombre de los archivos evolucionados me saque las graficas. Ya que son 5 archivos saca 5 graficas

def plots(file,i):
    
    #Importa los archivos producidos por $evolve.c$ y los separa
    input=np.loadtxt(file)
    t=input
    c1=t[:,1]
    c2=t[:,2]
    
    #Grafica los datos
    pylab.plot(c1,c2,'k*')
    pylab.xlabel('$Kpc$')
    pylab.ylabel('$Kpc$')
    pylab.title('Grafica posicion de Galaxias')
    
    #Guarda las graficas de los datos bajo el nombre del archivo importado con terminacion
    #.png
    alfa='Grafica'+str(i)
    pylab.savefig(alfa)
    pylab.close()

for i in range (len(sys.argv)):
    if i != 0:
        name = str(sys.argv[i])
        plots(name, i)