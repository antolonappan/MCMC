import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

rc('font', **{'family':'serif','serif':['Palatino'],'size':30})
rc('text', usetex=True)

d = np.loadtxt('chain1.csv', delimiter = ',')
d1 = np.loadtxt('chain2.csv', delimiter = ',')

chi = d[:,0]
a = d[:,1]
b = d[:,2]
c = d[:,3]

chi1 = d[:,0]
a1 = d1[:,1]
b1 = d1[:,2]
c1 = d1[:,3]
bins = np.linspace(0,1,15)
N=len(a)
N2=len(a1)
i = np.arange(0,N,1)
i1 = np.arange(0,N2,1)

#plt.hist(c1,bins)
#plt.plot(b)
#plt.plot(b1)
#plt.show()
#plt.scatter(i,b)
#plt.show()
#plt.scatter(i,c)
#plt.show()

af=a1[N/4:]
bf=b1[N/4:]
cf=c1[N/4:]
#df=d[N/2+N/4:]
#ef=e[N/2+N/4:]

ma = np.mean(af)
mb = np.mean(bf)
mc = np.mean(cf)
#md = np.mean(df)
#me = np.mean(ef)

sa = np.std(af) 
sb = np.std(bf)
sc = np.std(cf)
#sd = np.std(df)
#se=  np.std(ef)
k=b+c
k1=b1+c1
print(ma,mb,mc)
print(sa,sb,sc)
p=np.concatenate((a,a1),axis=1)
q=np.concatenate((k,k1),axis=1)
cov = np.cov(p,q )
lambda_, v = np.linalg.eig(cov)
lambda_ = np.sqrt(lambda_)
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
ax = plt.subplot(111)#, aspect='equal')
for j in xrange(1, 4):
    ell = Ellipse(xy=(np.mean(p), np.mean(q)),
                  width=lambda_[0]*j*2, height=lambda_[1]*j*2,
                  angle= np.rad2deg(np.arcsin(v[1, 0])))
    ell.set_facecolor('none')
    ax.add_artist(ell)
plt.scatter(p, q)
plt.xlabel(r'$H_o$')
plt.ylabel(r'$\Omega_K$')
plt.show()
