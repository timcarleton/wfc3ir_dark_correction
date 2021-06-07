from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import darkmodel2
import copy
import matplotlib.pyplot as plt

def getdarkcorrimg(original,refdir='/store/skysurf/',rawloc='find'):

    flt=fits.open(original)

    if 'N/A' not in flt[0].header['DARKFILE']:
        fdrk=fits.open(refdir+'crds_cache/references/hst/wfc3/'+flt[0].header['DARKFILE'].split('$')[-1])
        fdrkdata=fdrk['SCI'].data
    else:
        print('no dark file, correction not performed')
        return -1

        
    if 'N/A' not in flt[0].header['PFLTFILE']:
        fflt=fits.open(refdir+'crds_cache/references/hst/wfc3/'+flt[0].header['PFLTFILE'].split('$')[-1])
        ffltdata=fflt['SCI'].data
    else:
        ffltdata=np.ones_like(flt['SCI'].data)

        
    if 'N/A' not in flt[0].header['DFLTFILE']:
        dflt=fits.open(refdir+'crds_cache/references/hst/wfc3/'+flt[0].header['DFLTFILE'].split('$')[-1])
        dfltdata=dflt['SCI'].data
    else:
        dfltdata=np.ones_like(ffltdata)

    if rawloc=='find':
        yr=flt[0].header['DATE-OBS'].split('-')[0]
        rwloc='/media/WD10TB_cl1/skysurf/wfc3ir_rw/'+flt[0].header['FILTER']+'/'+yr+'/'
    else:
        rwloc=rawloc
        
    wd=np.where(fdrk['DQ'].data==0)

    root=original.split('/')[-1].split('_')[0]
    realdark=sigma_clipped_stats(fdrk['SCI'].data[wd])[0]*2.5/fdrk[0].header['EXPTIME']
    moddark=darkmodel2.darkmodelfile(original,rwloc+root+'_raw.fits')

    if flt['SCI'].data.shape[0]<1014:
        if 'N/A' not in flt[0].header['PFLTFILE']:
            offsety,offsetx=-int(flt[1].header['LTV1'])+5,-int(flt[1].header['LTV2'])+5
            ffltdat=ffltdata[offsetx:-offsetx,offsety:-offsety]
            dfltdat=dfltdata[offsetx:-offsetx,offsety:-offsety]
        else:
            ffltdat=ffltdata
            dfltdat=dfltdata
    else:
        if 'N/A' not in flt[0].header['PFLTFILE']:
            ffltdat=ffltdata[5:-5,5:-5]
            dfltdat=dfltdata[5:-5,5:-5]
        else:
            ffltdat=ffltdata
            dfltdat=dfltdata
        
    drk0=fdrk['SCI'].data[5:-5,5:-5]

    drkint=fits.open('/store/skysurf/reference/int_'+flt[0].header['SAMP_SEQ']+'.fits')
    drkm=fits.open('/store/skysurf/reference/slope_'+flt[0].header['SAMP_SEQ']+'.fits')

    drkdat=drkm[0].data*moddark+drkint[0].data

    newimg=flt['SCI'].data+(drk0*2.5/fdrk[0].header['EXPTIME']-drkdat)/ffltdat/dfltdat
    #newimg=flt['SCI'].data+2.5/ffltdat/dfltdat*drkdat/fdrk[0].header['EXPTIME']*(1-moddark/realdark)

    #newimg = flt['SCI'].data + (epsilon)*dark_frame # find epsilon to minimize quadrant differences
    newflt=copy.copy(flt)
    newflt['SCI'].data=newimg
    return newflt
