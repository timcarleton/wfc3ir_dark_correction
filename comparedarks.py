from astropy.io import fits,ascii
from scipy import stats
import numpy as np
from astropy.stats import sigma_clipped_stats
import darkmodel2
import matplotlib.pyplot as plt


plt.rcParams['xtick.top']=True
plt.rcParams['ytick.right']=True
plt.rcParams['xtick.direction']='in'
plt.rcParams['ytick.direction']='in'
plt.rcParams['font.family']='serif'
plt.rcParams['font.serif']=['Times']
plt.rcParams["figure.figsize"]=[16,12]
plt.rcParams['font.size']=30
plt.rcParams['image.cmap'] = 'plasma'
plt.rcParams['axes.linewidth']=4
plt.rcParams['axes.labelpad']=10
plt.rcParams['xtick.major.size']=9
plt.rcParams['ytick.major.size']=9
plt.rcParams['xtick.major.width']=3
plt.rcParams['ytick.major.width']=3
plt.rcParams['xtick.minor.size']=6
plt.rcParams['ytick.minor.size']=6
plt.rcParams['xtick.minor.width']=1
plt.rcParams['ytick.minor.width']=1
plt.rcParams['ytick.major.pad']=5
plt.rcParams['xtick.major.pad']=10
plt.rcParams['axes.labelsize']=35
plt.rcParams['axes.labelweight']='bold'
plt.rcParams['savefig.bbox']='tight'
plt.rcParams['savefig.dpi']=200
plt.rcParams['lines.linewidth']=4
plt.rcParams['lines.markersize']=10
plt.rcParams['errorbar.capsize']=0
plt.rcParams['hatch.linewidth']=2

dat=ascii.read('/store/skysurf/wfc3ir/F125W/wfc3ir_F125W_percentileclip_output_quadcorr.csv')
datgood=ascii.read('F125W_quad_imageinfo.csv')
print(datgood.keys())
#for i in range(len(dat)):

reals=[]
mods=[]
wgood=[]
delta=[]
delta2=[]
for i in range(len(dat)):
#for i in range(500):
    wi=np.where(dat['root'][i]==datgood['root'])[0]
    if len(wi)==0:
        continue
    wgood.append(i)
    flt=fits.open('/store/skysurf/wfc3ir/F125W/new_data/'+dat['root'][i]+'_new_flt.fits.gz')
    yr=flt[0].header['DATE-OBS'].split('-')[0]
    fdrk=fits.open('/home/tcn/crds_cache/references/hst/wfc3/'+flt[0].header['DARKFILE'].split('$')[-1])
    fflt=fits.open('/home/tcn/crds_cache/references/hst/wfc3/'+flt[0].header['PFLTFILE'].split('$')[-1])
    dflt=fits.open('/home/tcn/crds_cache/references/hst/wfc3/'+flt[0].header['DFLTFILE'].split('$')[-1])
    wd=np.where(fdrk['DQ'].data==0)
    reals.append(sigma_clipped_stats(fdrk['SCI'].data[wd])[0]*2.5/fdrk[0].header['EXPTIME'])
    moddark=darkmodel2.darkmodelfile('/store/skysurf/wfc3ir/F125W/new_data/'+dat['root'][i]+'_new_flt.fits.gz','/media/WD10TB_cl1/skysurf/wfc3ir_rw/F125W/'+yr+'/'+dat['root'][i]+'_raw.fits')
    mods.append(moddark)

    datleft=flt['SCI'].data[:,494:507]
    datright=flt['SCI'].data[:,507:520]

    wokl=np.where(flt['DQ'].data[:,494:507]==0)
    wokr=np.where(flt['DQ'].data[:,507:520]==0)
    
    delta.append(sigma_clipped_stats(datright[wokr])[0]-sigma_clipped_stats(datleft[wokl])[0])

    newimg=flt['SCI'].data-2.5/fflt['SCI'].data[5:-5,5:-5]/dflt['SCI'].data[5:-5,5:-5]*fdrk['SCI'].data[5:-5,5:-5]/fdrk[0].header['EXPTIME']*(1-moddark/reals[-1])
    #newimg=flt['SCI'].data-fdrk['SCI'].data[5:-5,5:-5]*(reals[-1]-moddark)

    newdatleft=newimg[:,494:507]
    newdatright=newimg[:,507:520]
    delta2.append(sigma_clipped_stats(newdatright[wokr])[0]-sigma_clipped_stats(newdatleft[wokl])[0])

wgood=np.array(wgood)
mods=np.array(mods)
reals=np.array(reals)
delta=np.array(delta)
delta2=np.array(delta2)


plt.clf()
#plt.plot(reals-mods,abs(dat['q1_corr'][wgood])+abs(dat['q2_corr'][wgood])+abs(dat['q3_corr'][wgood])+abs(dat['q4_corr'][wgood]),'o')
plt.plot(reals-mods,delta,'o',alpha=.25,label='No Correction')
plt.plot(reals-mods,delta2,'o',label='With Dark Adjustment')

plt.ylim(-.01,.01)
plt.legend()
plt.xlabel('Calibration Dark Level - Model Dark Level (e/s)')
plt.ylabel('left border - right border (e/s)')
plt.savefig('quadvsdiff.png')

plt.clf()

