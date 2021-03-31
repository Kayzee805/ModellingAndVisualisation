import cythonCode 

n=100
epsilon=0.001
method="gauss"
#cythonCode.runMagnetic(n,method,epsilon)
cythonCode.electricField(n,method,epsilon)