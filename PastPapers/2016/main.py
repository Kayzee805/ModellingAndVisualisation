import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
from IPython.core.display import clear_output
from tqdm import tqdm
from scipy.stats import sem


class Model(object):

    def __init__(self):
        super().__init__()