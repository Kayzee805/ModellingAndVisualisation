import numpy as np
import time
np.random.seed(10)


x=np.zeros((50,50,50))
noise = np.random.uniform(-1,1,(50,50,50))
# x=np.add(x,noise)
test=(np.indices((50,50,50)).sum(axis=0)%2)

x[test==1]=5
a = x[test==1]
print(a.shape)
