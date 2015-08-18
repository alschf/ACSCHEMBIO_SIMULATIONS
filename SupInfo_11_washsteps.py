import numpy

def getChangeComplex(concLt, concP, n, secs, koff, kon):
    concL = 0.0
    for s in range(0, secs):
        for i in range(0,n):
            dt = 1/float(n)
            dPL = -koff*(concLt - concL)*dt + kon*concP*concL*dt
            concL = concL - dPL
            ##if i%100000 == 0:
            ##    print "dt is" , dt
            ##    print dPL
            ##    print concPL
            ##    print concL
    return (concLt-concL)/concLt
    
    
getChangeComplex(1E-15, 1E-6, 100000, 1, 1E0, float(1E6))
