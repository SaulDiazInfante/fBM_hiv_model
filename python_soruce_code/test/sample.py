# flake8: noqa
from fbm import FBM
import matplotlib.pyplot as plt
import time
import math


def h(t):
    # return 0.499*math.sin(t) + 0.5
    # return 0.6 * t + 0.3
    return 0.5 * math.exp(-8*t**2) + 0.35


fbm_generator = FBM(2 ** 8, .88)
t = fbm_generator.times()
fbm_realization = fbm_generator.fbm()
fgn_realization = fbm_generator.fgn()
plt.plot(t, fbm_realization)
plt.plot(t[0:-1], fgn_realization)
plt.show()
