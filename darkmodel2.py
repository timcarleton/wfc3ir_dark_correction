from astropy.io import fits
import numpy as np
import getoverscanslope
import geti087volt

def getfitparams(sampseq):
    if sampseq=='SPARS200':
        fit=[5.28690098, -0.30673704, -5.19573127]
    elif sampseq=='SPARS100':
        fit=[6.57510813,  2.73692116, -2.46815561]
    elif sampseq=='SPARS50':
        fit=[7.09722715,  2.54281641, -2.33913029]
    elif sampseq=='SPARS25':
        fit=[4.99023146,  0.12734818, -4.34830231]
    elif sampseq=='SPARS10':
        fit=[3.19570852,  0.05296631, -1.14721234]
    elif sampseq=='RAPID':
        fit=[1.53886565e+00, 1.32664978e-03, 4.88868024e-01]
    elif sampseq=='STEP25':
        fit=[3.23120847,  0.21606371, -2.18638573]
    elif sampseq=='STEP50':
        fit=[6.96449495,  5.87656204, -0.40520036]
    elif sampseq=='STEP100':
        fit=[6.25962986,  1.99431816, -2.63637995]
    elif sampseq=='STEP200':
        fit=[6.36874403,  2.19981356, -2.94320471]
    elif sampseq=='STEP400':
        fit=[7.02796073, 5.79328485, 1.42697463]
    return fit

def darkmodelfile(fltfile,rawfile):

    slope=getoverscanslope.getoverscanslope(rawfile)
    f=fits.open(fltfile)
    
    voltsi=geti087volt.getvolt(f[0].header['EXPSTART']+np.linspace(f[0].header['SAMPZERO'],f[0].header['EXPTIME'],f[0].header['NSAMP']-1)/60/60/24)
    
    deltavolt=np.polyfit(np.linspace(f[0].header['SAMPZERO'],f[0].header['EXPTIME'],f[0].header['NSAMP']-1),voltsi,1)[0]

    fit=getfitparams(f[0].header['SAMP_SEQ'])
    
    return fit[0]/100+fit[1]*slope+fit[2]*deltavolt

def darkmodel(overscanslope,deltavolt,sampseq):
    if type(sampseq)==str:
        fit=getfitparams(sampseq)
        return fit[0]/100+fit[1]*overscanslope+fit[2]*deltavolt
    else:
        res=np.zeros(len(sampseq))+np.nan

        exsamps=['SPARS200','SPARS100','SPARS50','SPARS25','SPARS10','RAPID','STEP400','STEP200','STEP100','STEP50','STEP25']
        for i in exsamps:
            wi=np.where(sampseq==i)[0]
            fiti=getfitparams(i)
            res[wi]=fiti[0]/100+fiti[1]*overscanslope[wi]+fiti[2]*deltavolt[wi]

        return res


def getparam(fltfile,rawfile):
    slope=getoverscanslope.getoverscanslope(rawfile)
    deltavolt=[]
    f=fits.open(fltfile)
    
    voltsi=geti087volt.getvolt(f[0].header['EXPSTART']+np.linspace(f[0].header['SAMPZERO'],f[0].header['EXPTIME'],f[0].header['NSAMP']-1)/60/60/24)
    deltavolt=np.polyfit(np.linspace(f[0].header['SAMPZERO'],f[0].header['EXPTIME'],f[0].header['NSAMP']-1),voltsi,1)[0]

    return slope,deltavolt
