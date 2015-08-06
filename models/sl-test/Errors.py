#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

#dt   = [0.015708,0.00785398,0.00392699,0.0019635]
#err1 = [0.0338165,0.0342615,0.0344928,0.0345968]
#err8 = [0.000267513,8.69884e-05,6.3132e-05,6.59601e-05]

#dt  = [0.00785398,0.00392699,0.0019635]
#err = [0.00849903,0.00943511,0.00981663]

#dx  = [2.0/64,2.0/128,2.0/256]
#err = [0.0616072,0.00943511,0.0011498]
dx  = [2.0/64,2.0/128,2.0/256,2.0/512]
err = [0.06221,0.0096965,0.00123947,0.000166176]

plt.loglog(dx,err,'o-')
plt.show()

log_x = np.log(dx)
log_e = np.log(err)

y   = np.vstack([log_x,np.ones(len(log_x))]).T
m,c = np.linalg.lstsq(y,log_e)[0]
print m

