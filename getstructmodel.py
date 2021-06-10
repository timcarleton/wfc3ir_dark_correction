import numpy as np
from astropy.io import fits
import darkmodel2


def getstructmodel(fltfile,rawfile,sptfile='',param4=False):

    f=fits.open(fltfile)
    slope,deltavolt=darkmodel2.getparam(fltfile,rawfile)
    
    param0file=fits.open('/store/skysurf/reference/param0_'+f[0].header['SAMP_SEQ']+'.fits')
    param1file=fits.open('/store/skysurf/reference/param1_'+f[0].header['SAMP_SEQ']+'.fits')
    param2file=fits.open('/store/skysurf/reference/param2_'+f[0].header['SAMP_SEQ']+'.fits')
    
    return param0file[0].data/100+param1file[0].data*slope+param2file[0].data*deltavolt
