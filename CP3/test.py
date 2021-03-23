import numpy as np
np.random.seed(5)
x = np.ones((3,3,3))
noise = np.random.uniform(-1,1,(3,3,3))
x=np.add(x,noise)

xgrad = np.gradient(x)

# print(x)
# print("hello")

# mid = 0
# manual = x[mid,mid,(mid+1)%3]+x[mid,mid,(mid-1)%3]+ x[(mid+1)%3,mid,mid]+x[(mid-1)%3,mid,mid]+x[mid,(mid+1)%3,mid]+x[mid,(mid-1)%3,mid]

# roll = np.roll(x,1,axis=0)+np.roll(x,-1,axis=0)+np.roll(x,1,axis=1)+np.roll(x,-1,axis=1)+np.roll(x,1,axis=2)+np.roll(x,-1,axis=2)
# auto = roll[mid,mid,mid]

# print(f"Manual={manual}  auto={auto}")


print(0.111111111111111111119)
x=np.array(0.111111111111111111119)
print(x)
# print(f"rolling in x={np.roll(x,1,axis=0)}")
# print(f"rolling in y={np.roll(x,1,axis=1)}")
# print(f"rolling in z={np.roll(x,1,axis=2)}")