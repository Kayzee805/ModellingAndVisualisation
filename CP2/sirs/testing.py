import logging
import threading
import time
import numpy as np
from astropy.stats import jackknife_resampling
from astropy.stats import jackknife_stats



x = np.linspace(0.2,0.5,16)
print(x)