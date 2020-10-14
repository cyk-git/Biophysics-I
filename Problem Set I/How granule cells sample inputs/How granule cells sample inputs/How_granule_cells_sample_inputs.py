import numpy as np 
from matplotlib import pyplot as plt 
 
M=15
N=45
KMax = M
K = np.arange(0,KMax+1) 
p = np.zeros([KMax+1],'f4')
for n in range(KMax):
    p[n] = 1
    for i in range(N):
        c = 1
        for j in range(K[n]):
            c*=((M-j)/(j+1))
        p[n]*=(1-i/c)

plt.title("Probability") 
plt.xlabel("K") 
plt.ylabel("p") 
plt.plot(K,p) 
plt.show()
