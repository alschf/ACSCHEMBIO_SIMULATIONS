import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy as sy
import pylab as plb  

## for ka = 6.45 and yield of 45%
x = numpy.array([5*10**-8, 5*10**-7, 5*10**-6])
y1= [6.21211548899e-14,6.47203509909e-12, 2.59504681254e-11]


## for ka = 6.45 and yield of 80%
y2= [9.02760626534e-14,1.04469919517e-11,4.61341655564e-11]

## for ka = 6.1 and yield of 80%
y3= [1.14736906745e-14,3.21005820379e-12,3.61926641064e-11]

## for ka = 6.1 and yield of 45%
y4= [6.45395100445e-15,1.80565773963e-12,2.03583735598e-11]

ylist = [y1, y2, y3, y4]

def func(x, Y, k):
  '''fits EQ4'''
  del_conc = 5*10**-5
  del_diversity = 707788
  return (((del_conc)/del_diversity)*Y*( k**3) * (x**3))/(1 + 3*k*x + 3* (k**2) * (x**2) + (k**3)*(x**3))

##p0 provides an estimated yield and Ka value
p0 = sy.array(([0.2,7*10**6]))
coefflist = []

for u in ylist:
    coeffs, matcov = curve_fit(func, x, u, p0)
    coefflist.append(coeffs)


xplot = np.linspace(5*10**-9, 5*10**-6, 1000)

fitlist = []
for u in coefflist:
    yaj = func(xplot, u[0], u[1])
    fitlist.append(yaj)
    print(u)


plt.plot(x,y1,'x',x,y2,'o',xplot,fitlist[0],'r-',xplot,fitlist[1],'b-',xplot,fitlist[2],'y-',xplot,fitlist[3],'g-')
plt.xlabel('[Protein]')
plt.ylabel('final [Ligand]')
plt.show()
