# Installation

python setup.py install --user

# usage: example for a 2 layer problem

```
import def_radius
import numpy as np

H = np.array([500, 3000]) # Layer thickness (m)
gp = np.array([1e-2])     # reduced gravity (m/s^2) -- density jumps between layers
f = 1e-4                  # coriolis parameter (1/s)

rd = def_radius.cal_rad(H,gp,f)

print("first deformation radius: {0:.1f} km".format(rd[1]*1e-3))
```