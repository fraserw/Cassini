#! /usr/bin/env python

import numpy as num, pylab as pyl,sys
from scipy import interpolate as interp
import scipy as sci
from numpy import linalg
from astropy.io import fits as pyf

viswave=num.arange(96)*0.0073+0.35+0.00036
irwave=num.arange(256)*0.0166+0.85+0.0083
def getSpec(imData,l,s,cosmoCut=1.8,channelCut=256):
    w=num.arange(256)*0.0166+0.85+0.0083
    #trim known bad spaxels
    badChannels=[22,45,46,47,48,97,122,180,205,208,211,226]

    spec=imData[:,l,s]*1.
    for jj in range(len(badChannels)):
        #spec[badChannels[jj]-1:badChannels[jj]+2]=num.nan
        spec[badChannels[jj]]=num.nan
    #now look for cosmic rays
    for jj in range(len(spec)):
        ww=num.where(num.isnan(spec[max(0,jj-3):min(255,jj+4)])<>True)
        
        m=num.median(spec[max(0,jj-3):min(255,jj+4)][ww])
        std=num.std(spec[max(0,jj-3):min(255,jj+4)][ww])
        #print jj,spec[jj],m,std,spec[max(0,jj-3):min(255,jj+4)][ww]
        if abs(spec[jj]-m)/std>cosmoCut:
            #print jj,spec[jj],m
            spec[jj]=m#num.nan
    return (w[:channelCut],spec[:channelCut])

def getSpecVis(imData,l,s,cosmoCut=1.75,channelCut=96):
    w=num.arange(96)*0.0073+0.35+0.00036
    #trim known bad spaxels
    badChannels=[]

    spec=imData[:,l,s]*1.
    for jj in range(len(badChannels)):
        #spec[badChannels[jj]-1:badChannels[jj]+2]=num.nan
        spec[badChannels[jj]]=num.nan
    #now look for cosmic rays
    for jj in range(len(spec)):
        ww=num.where(num.isnan(spec[max(0,jj-3):min(95,jj+4)])<>True)
        
        m=num.median(spec[max(0,jj-3):min(95,jj+4)][ww])
        std=num.std(spec[max(0,jj-3):min(95,jj+4)][ww])
        
        #print jj,spec[jj],m,std,spec[max(0,jj-3):min(95,jj+4)][ww]
        #if std==0:
        #    print spec
        #    sys.exit()
        if abs(spec[jj]-m)/std>cosmoCut:
            #print jj,spec[jj],m
            spec[jj]=m#num.nan
    return (w[:channelCut],spec[:channelCut])


def water(wave,spec,returnSpec=False):
    w1=num.where((wave>=1.28)&(wave<=1.40))
    w2=num.where((wave>=1.7)&(wave<=1.85))
    w3=num.where((wave>=2.2)&(wave<=2.26))
    w4=num.where((wave>=3.52)&(wave<=3.75))

    medWaves=num.array([num.nanmedian(wave[w1]),num.nanmedian(wave[w2]),num.nanmedian(wave[w3]),num.nanmedian(wave[w4])])
    medSpec=num.array([num.nanmedian(spec[w1]),num.nanmedian(spec[w2]),num.nanmedian(spec[w3]),num.nanmedian(spec[w4])])

    f=interp.interp1d(medWaves,medSpec)

    wf1=num.where((wave>=1.48)&(wave<=1.55))
    wf2=num.where((wave>=1.95)&(wave<=2.05))
    wf3=num.where((wave>=2.8)&(wave<=3.03))

    c1=f(num.nanmedian(wave[wf1]))
    c2=f(num.nanmedian(wave[wf2]))
    c3=f(num.nanmedian(wave[wf3]))

    d1=num.max((c1-num.nanmedian(spec[wf1]))/c1,0.0)
    d2=num.max((c2-num.nanmedian(spec[wf2]))/c2,0.0)
    d3=num.max((c3-num.nanmedian(spec[wf3]))/c3,0.0)
    
    #w=num.where((wave>=medWaves[0])&(wave<=medWaves[2]))
    #featureSpec=spec[w]-f(wave[w])
    #featureWave=wave[w]
    #pyl.plot(wave[w],featureSpec)
    #pyl.title(d1+d2)
    #pyl.show()
    if returnSpec:
        www=num.where((wave>num.min(medWaves))&(wave<num.max(medWaves)))
        return(wave[www],spec[www]-f(wave[www]),d1+d2,d3)
    return (d1+d2,d3)

def oSlope(wave,spec):
    w=num.where((wave>=0.5)&(wave<=0.7))

    A=[]
    for ii in range(len(w[0])):
        A.append([1.,wave[w[0][ii]]])
    A=num.array(A)
    At=A.transpose()
    X=num.dot(num.dot(linalg.inv( num.dot(At,A) ), At), spec[w])

    l=X[0]+X[1]*wave[w]
    reddening=(l[len(l)-1]-l[0])/l[0]/2.
    return reddening
    #pyl.plot(wave[w],l)
    #pyl.plot(wave,spec)
    #pyl.show()
    #sys.exit()

def dSlope(wave,spec):
    w=num.where((wave>=0.36)&(wave<=0.48))

    A=[]
    for ii in range(len(w[0])):
        A.append([1.,wave[w[0][ii]]])
    A=num.array(A)
    At=A.transpose()
    X=num.dot(num.dot(linalg.inv( num.dot(At,A) ), At), spec[w])

    l=X[0]+X[1]*wave[w]
    reddening=(l[len(l)-1]-l[0])/l[0]/1.2
    return reddening
    #pyl.plot(wave[w],l)
    #pyl.plot(wave,spec)
    #pyl.show()
    #sys.exit()

def silicate(wave,spec,returnSpec=False):
    w1=num.where((wave>=0.85)&(wave<=1.15))
    w2=num.where((wave>=2.18)&(wave<=2.3))

    m1=num.nanmedian(spec[w1])
    m2=num.nanmedian(spec[w2])
    return (m2-m1)/m2

    
if __name__=="__main__":
    import sys


    with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/2004163T121836_2004163T192848/cv1465662631_1_vis_strpd.fits') as han:
        data=han[0].data
    (channels,A,B)=data.shape
    for l in range(A):
        for s in range(B):
            (w,spec)=getSpecVis(data,l=l,s=s,channelCut=channels-1,cosmoCut=1.75)
            if num.median(spec)>0:
                pyl.plot(w,spec)
        #pyl.title(l)
    pyl.show()
    sys.exit()

    with pyf.open('cv1465662631_1_ir.fits') as han:
        data=han[0].data
    (channels,A,B)=data.shape

    n=0
    waterDepth=[]
    for l in range(A):
        waterDepth.append([])
        for s in range(B):
            (w,spec)=getSpec(data,l=l,s=s)
            if num.nanmedian(spec)>0.005:
                (contW,contSpec,wd)=water(w,spec,returnSpec=True)
                sd=silicate(w,spec)
                #waterDepth[len(waterDepth)-1].append(water(w,spec))
                #pyl.plot(w,spec)
                pyl.plot(w,spec)
                n+=1
                #if sd>0.3:
                #    pyl.plot(w,spec)
                #if wd>0.4:
                    #pyl.plot(contW,contSpec,'k-',label=str(l)+' '+str(s)+' '+str(wd))
                #if wd<0.1:
                    #pyl.plot(contW,contSpec,'r-',label=str(l)+' '+str(s)+' '+str(wd))
            else:
                waterDepth[len(waterDepth)-1].append(0.0)
        #if n>50:break
    #pyl.legend()
    pyl.show()
    sys.exit()
    pyf.writeto('junk.fits',num.array(waterDepth),clobber=True)
