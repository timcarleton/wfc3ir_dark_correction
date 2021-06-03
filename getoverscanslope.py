import numpy as np
from astropy.io import fits

def getoverscanslope(rawfile):
    f=fits.open(rawfile)
    means=[]
    times=[]
    
    for i in range(1,f[0].header['NSAMP']):
        mean1=np.mean(f['SCI',i].data[1:5,5:-5])
        mean2=np.mean(f['SCI',i].data[-5:-1,5:-5])
        means.append(np.mean([mean1,mean2]))
        times.append(f['SCI',i].header['SAMPTIME'])

    fit=np.polyfit(times,means,1)
    return fit[0]
