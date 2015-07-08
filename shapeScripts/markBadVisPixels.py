#! /usr/bin/env python

import numdisplay
from astropy.io import fits as pyf
import numpy as num,pylab as pyl
import sys,os
from astLib import astImages as Im

def button(event):
    global bads,badLines
    x,y=event.xdata,event.ydata
    if event.button==1:
        bads.append([int(x+0.5),int(y+0.5)])
    elif event.button==3:
        badLines.append(int(x+0.5))
        
def getScaledImage(sec,A=501,B=501,nsamp=200,contrast=0.4,lims=[]):
    if lims==[]:
        lims=[0,A,0,B]
    im=num.zeros([A,B]).astype(float)
    (z1,z2)=numdisplay.zscale.zscale(sec,nsamples=nsamp,contrast=contrast)
    norm=Im.normalise(sec,[z1,z2])
    #(a,b)=sec.shape                                                            
    #im[:a,:b]=norm                                                             
    im[lims[0]:lims[1],lims[2]:lims[3]]=norm
    return im

#input like cv1465669068_1_vis_mean.fits

if len(sys.argv)==1:
    print 'Provide an image name like cv1465669068_1_vis_mean.fits'
    sys.exit()
imName=sys.argv[1]
mapName=imName.split('.')[0]+'_bpmap.fits'

if len(sys.argv)>2:
    contrast=float(sys.argv[2])
else:
    contrast=2.5
bads=[]
badLines=[]

with pyf.open(imName) as han:
    data=han[0].data

(A,B)=data.shape
print A,B
scaled=getScaledImage(data,A=A,B=B,contrast=contrast)

fig=pyl.figure()
sp1=fig.add_subplot(1,1,1)
sp1.set_aspect('equal','datalim')
sp1.autoscale_view(True,True,True)

axIm=pyl.imshow(scaled*1.,cmap='gray',interpolation='nearest')

pyl.connect('button_press_event',button)
pyl.show()

for i in range(len(badLines)):
    for j in range(A):
        bads.append([badLines[i],j])
bads=num.array(bads)

for i in range(len(bads)):
    if bads[i][1]<A and bads[i][0]<B:
        data[bads[i][1],bads[i][0]]=0.0

newscaled=getScaledImage(data,A=A,B=B,contrast=2.5)

fig=pyl.figure()
sp1=fig.add_subplot(1,1,1)
sp1.set_aspect('equal','datalim')
sp1.autoscale_view(True,True,True)

axIm=pyl.imshow(newscaled*1.,cmap='gray',interpolation='nearest')

pyl.connect('button_press_event',button)
pyl.show()

w=num.where(data<>0)
data[w]=1.

HDU=pyf.PrimaryHDU(data)
List=pyf.HDUList([HDU])
try:
    os.remove(mapName)
except: pass
List.writeto(mapName)

