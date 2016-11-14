#!/usr/bin/env python3
import numpy as np
cubicfit = np.array([-0.0625, 0.375, -1.25, 0.625, 0.3125])
ls3rdorder = np.array([0, 0.166666667, -1, 0.5, 0.3333333333])
ls4thorder = np.array([-0.083, 0.5, -1.5, 0.8333333333, 0.25])

print(ls3rdorder - cubicfit)
print(ls4thorder - cubicfit)

