import numpy as np
dat=np.loadtxt('/store/skysurf/reference/I087.txt')

def getvolt(jd):
    i=np.searchsorted(dat[:,0],jd)
    return dat[i-1,1]
    
def getdeltavolt(jdstart,jdend):
    i=np.searchsorted(dat[:,0],jdstart)
    j=np.searchsorted(dat[:,0],jdend)
    fit=np.polyfit(dat[i-1:j-1,0],dat[i-1:j-1,1],1)
    return fit[0]
