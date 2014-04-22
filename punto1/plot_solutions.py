#codigo que grafica los datos guardados en datos xy, archivo creado anteriormente

import numpy as np, pylab

datos = np.loadtxt('datosxy.dat')

for n in range(30):
    x0 =datos[999*n:998+999*n,0]
    y0 =datos[999*n:998+999*n,1]
    pylab.plot(x0,y0,'-b')
    pylab.xlabel('$x$')
    pylab.ylabel('$y$')
    pylab.title('$x_0 = '+str(x0[0])+'\ \ y_0 = '+str(y0[0])+'$')
    pylab.savefig('volterra_'+str(30-n)+'.png',dpi=200)
    pylab.close()
