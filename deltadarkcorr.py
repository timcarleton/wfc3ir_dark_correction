from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
import darkmodel2
import copy

def getdarkcorrimg(original,refdir='/store/skysurf/',rawloc='find'):

    flt=fits.open(original)

    fdrk=fits.open(refdir+'crds_cache/references/hst/wfc3/'+flt[0].header['DARKFILE'].split('$')[-1])
    fflt=fits.open(refdir+'crds_cache/references/hst/wfc3/'+flt[0].header['PFLTFILE'].split('$')[-1])
    dflt=fits.open(refdir+'crds_cache/references/hst/wfc3/'+flt[0].header['DFLTFILE'].split('$')[-1])

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
        offsety,offsetx=-int(flt[1].header['LTV1'])+5,-int(flt[1].header['LTV2'])+5
        ffltdat=fflt['SCI'].data[offsetx:-offsetx,offsety:-offsety]
        dfltdat=dflt['SCI'].data[offsetx:-offsetx,offsety:-offsety]
    else:
        ffltdat=fflt['SCI'].data[5:-5,5:-5]
        dfltdat=dflt['SCI'].data[5:-5,5:-5]
        
    drkdat=fdrk['SCI'].data[5:-5,5:-5]

    newimg=flt['SCI'].data-2.5/ffltdat/dfltdat*drkdat/fdrk[0].header['EXPTIME']*(1-moddark/realdark)

    newflt=copy.copy(flt)
    newflt['SCI'].data=newimg
    return newflt
