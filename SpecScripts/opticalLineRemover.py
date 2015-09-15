#! /usr/bin/env python

import astropy.io.fits as pyf,numpy as num,pylab as pyl
import optparse,numdisplay
from astLib import astImages as Im

def getEmptyRows(event):
    global emptyRows
    x,y=int(event.xdata+0.5),int(event.ydata+0.5)
    emptyRows.append(y)
    return
def getEmptyColumns(event):
    global emptyColumns
    x,y=int(event.xdata+0.5),int(event.ydata+0.5)
    emptyColumns.append(x)
    return
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

parser=optparse.OptionParser()
(opt,args)=parser.parse_args()

 
imageName=args[0]#'2004163T121836_2004163T192848/cv1465665036_1'
with pyf.open('/data/VIMS/covims_0004/procdata/'+imageName+'_vis.fits',ignore_missing_end=True) as han:
    imDat=han[0].data
    header=han[0].header

print imDat.shape
with pyf.open('/data/VIMS/covims_0004/procdata/'+imageName+'_vis_mean.fits') as han:
    meanDat=han[0].data

with open('/data/VIMS/covims_0004/procdata/'+imageName+'.parFile') as han:
    data=han.readlines()
    visCut=float(data[0].split()[0])


fig=pyl.figure(1)
sp1=fig.add_subplot(111)
pyl.imshow(num.clip(meanDat,max(num.min(meanDat),0.0),visCut),cmap='gray',origin='lower',interpolation='nearest')
pyl.title('Zoom to the region containing Phoebe.')
#pyl.connect('button_press_event',getEmptyRows)
pyl.show()

(x0,x1)=sp1.get_xlim()
(y0,y1)=sp1.get_ylim()


(N,A,B)=imDat.shape
emptyRows,emptyColumns=[],[]
for i in range(A):
    if i<y0 or i>y1:
        emptyRows.append(i)
for i in range(B):
    if i<x0 or i>x1:
        emptyColumns.append(i)

if len(emptyRows)<>0:
    for j in range(N):
        for i in range(B):
            meanStrip=num.median(imDat[j,:,i][emptyRows])
            imDat[j,:,i]-=meanStrip
    for i in range(B):
        meanStrip=num.median(meanDat[:,i][emptyRows])
        meanDat[:,i]-=meanStrip

#emptyColumns=[]
#scaled=getScaledImage(meanDat,A=A,B=B,contrast=2.)
#pyl.imshow(scaled,cmap='gray',origin='lower',interpolation='nearest')
#pyl.title('Select rows left and right of the target that do not contain Phoebe flux.')
#pyl.connect('button_press_event',getEmptyColumns)
#pyl.show()
#
#emptyColumns=num.unique(num.array(emptyColumns))
for j in range(N):
    for i in range(A):

        meanStrip=num.median(imDat[j,i,:][emptyColumns])
        imDat[j,i,:]-=meanStrip
for i in range(A):
    meanStrip=num.median(meanDat[i,:][emptyColumns])
    meanDat[i,:]-=meanStrip


pyf.writeto('/data/VIMS/covims_0004/procdata/'+imageName+'_vis_strpd.fits',imDat,header=header,clobber=True)
pyf.writeto('/data/VIMS/covims_0004/procdata/'+imageName+'_vis_mean_strpd.fits',meanDat,header=header,clobber=True)
    
        
#70650
#79413
    
