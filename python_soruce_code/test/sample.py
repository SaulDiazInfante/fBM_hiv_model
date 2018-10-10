# flake8: noqa
from fbm import FBM
import matplotlib.pyplot as plt
import time
import math
import numpy as np


def h(s):
    # return 0.499*math.sin(t) + 0.5
    # return 0.6 * t + 0.3
    return 0.5 * np.exp(-8.0 * s ** 2)


fbm_generator = FBM(2 ** 8, 0.85)
t = fbm_generator.times()
fbm_realization = fbm_generator.fbm()
fgn_realization = fbm_generator.fgn()
h_t = np.array(h(fbm_realization))
plt.plot(t, fbm_realization)
plt.plot(t, h_t)
# plt.plot(t[0:-1], fgn_realization)
plt.show()
