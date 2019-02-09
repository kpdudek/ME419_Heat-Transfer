#!/usr/bin/env

import numpy as np
import matplotlib.pyplot as mp
from matplotlib import rc

#rc('text',usetex=True)


###- PLOT 1 -###
Lf = np.linspace(0,.001,num=15)

q_dot = ((60.0-20.0)/((1.0/50.0)+(Lf/.025))) + ((60.0-30.0)/(.001/.05))

#print(q_dot)

mp.figure()
mp.plot(Lf,q_dot)
mp.xlabel('Film Thickness (m)')
mp.ylabel(r'Heat Flux $(\frac{W}{m^2})$')
mp.title(r'Transparent Film Thickness vs Heat Flux')



###- PLOT 2 -###
q_dot2 = (((60000*Lf+60)-20.0)/((1.0/50.0))) + (((60000*Lf+60)-30.0)/((Lf/.025)+(.001/.05)))


mp.figure()
mp.plot(Lf,q_dot2)
mp.xlabel('Film Thickness (m)')
mp.ylabel(r'Heat Flux $(\frac{W}{m^2})$')
mp.title(r'Opaque Film Thickness vs Heat Flux')


# Display figure windows
mp.show()

































