import numpy as np

x = np.ones((5,5,5))
noise = np.random.uniform(-1,1,(5,5,5))
x=np.add(x,noise)

xgrad = np.gradient(x)

print(np.shape(xgrad))
print(xgrad)