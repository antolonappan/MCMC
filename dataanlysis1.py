import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.patches import Ellipse
from matplotlib import rc

rc('font', **{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)

dat=np.loadtxt('m.csv',delimiter =',')

X = dat[:,0]
Y = dat[:,1]
dy = dat[:,2]
plt.errorbar(X,Y,xerr=0, yerr=dy,fmt='o',label='Data')
H=[]

for i in X:
     H = np.append(H,np.sqrt(pow(70,2)*(.3*pow((1+i),3)+.7)) )  
     

plt.scatter(X,H,label='Fitted',c='red')
plt.xlabel(r'$z$')
plt.ylabel(r'$H(z)$')
legend = plt.legend(loc='upper center', shadow=True)
plt.show()


d = np.loadtxt('chain1.csv', delimiter = ',')

chi = d[:,0]
a = d[:,1]
b = d[:,2]
c = d[:,3]
#d = d[:,4]
#d1 = np.loadtxt('chain1.csv', delimiter = ',')
#e = d1[:,5]

N=len(a)
#i = np.arange(0,N,1)
#plt.scatter(i,a)
#plt.show()
#plt.scatter(i,b)
#plt.show()
#plt.scatter(i,c)
#plt.show()

af=a[N/2+N/4+N/8:]
bf=b[N/2+N/4+N/8:]
cf=c[N/2+N/4+N/8:]
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

print(ma,mb,mc)
print(sa,sb,sc)

y1=[]
for x in X:
       y1 = np.append(y1,(ma*pow(x,2))+(mb*x)+mc)


#plt.plot(X,y)
#plt.show()

#rc('font', **{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)


f, axarr = plt.subplots(2, 2)

axarr[0, 0].scatter(a, chi,s=.5,c='r')
axarr[0, 0].set_title('a')

axarr[0, 1].scatter(b, chi,s=.5,c='r')
axarr[0, 1].set_title('b')
axarr[1, 0].scatter(c, chi,s=.5,c='r')
axarr[1, 0].set_title('c')
#axarr[1, 1].scatter(d, chi,s=.5,c='r')
#axarr[1, 1].set_title('chi vs d')

plt.figure()
plt.errorbar(X,Y,xerr=0, yerr=dy)

plt.plot(X,y1)
plt.show()

cov = np.cov(a, b)
lambda_, v = np.linalg.eig(cov)
lambda_ = np.sqrt(lambda_)
from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
ax = plt.subplot(111, aspect='equal')
for j in xrange(1, 4):
    ell = Ellipse(xy=(np.mean(x), np.mean(y)),
                  width=lambda_[0]*j*2, height=lambda_[1]*j*2,
                  angle=np.rad2deg(np.arccos(v[0, 0])))
    ell.set_facecolor('none')
    ax.add_artist(ell)
plt.scatter(x, y)
plt.show()

