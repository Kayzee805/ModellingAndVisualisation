import numpy as np
import random
import time 

t1=time.time()
for i in range(1000000):
    x = random.randint(0,51)
t2=time.time()
for i in range(1000000):
    y = np.random.randint(0,50)
t3=time.time()

print(f"Random = {t2-t1} numpy = {t3-t2}")