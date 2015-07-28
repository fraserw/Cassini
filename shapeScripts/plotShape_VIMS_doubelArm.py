#! /usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyl
from matplotlib import backend_bases
import numpy as num, scipy as sci
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys,os
import illuminate as ill
import pixelLatLong as pll
import pickle as pick
import random
from astropy.io import fits as pyf
import time
import emcee
from scipy import signal
from pixTimes import getPixTimes,getVels
from MChandler import getFit
from shapeLib import *
from lineFitUtils import poly
import specAnalysis
#from numpy import linalg

def getPars(fn):
    with open(fn) as han:
        r=han.readlines()
    vc=float(r[0].split()[0])
    vi=float(r[1].split()[0])
    s=r[2].split()[0].split(',')
    oxv=float(s[0])
    oyv=float(s[1])
    s=r[3].split()[0].split(',')
    oxi=float(s[0])
    oyi=float(s[1])
    try:
        s=r[4].split()
        dn=float(s[0])
    except:
        dn=None
    return (vc,vi,oxv,oyv,oxi,oyi,dn)


def callShapeGen(r,vertices,vertIndices,ps_vis,ps_ir,inot,jnot,dXV,dYV,dZV,dXI,dYI,dZI,Ls,ls,pTimesV,pTimesI,imData,az_adjust=0.0):
    global steps
    steps+=1
    print steps,

    (L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi,ov,ova)=r

    if abs(oxv)>1200. or abs(oyv)>1200. or oxv==num.inf or oxv==-num.inf or oyv==num.inf or oyv==-num.inf:
        #print oxv,oyv
        return -num.inf
    if abs(oxi)>1200. or abs(oyi)>1200. or oxi==num.inf or oxi==-num.inf or oyi==num.inf or oyi==-num.inf:
        #print oxi,oyi
        return -num.inf
    if  dn<0 or ov<0 or ov>5 or ova<0 or ova>=360. or dn<0:
        #print dn,ov,ova
        return -num.inf
    if L_o_n>=360. or L_o_n<0. or l_o_n>90. or l_o_n<-90. or A_o_n>=360. or A_o_n<0. or l_s>90. or l_s<-90. or L_s>=360. or L_s<0. or A_s<0 or A_s>=360.:
        print L_s,l_s,A_s
        return -num.inf

    print L_o_n,l_o_n,A_o_n
    print '    ',L_s,l_s,A_s
    print '    ',dn
    print '    ',oxv,oyv,oxi,oyi
    print '    ',ov,ova
    (imageVis,poly3d,colours,rot_vert,vertsInPixelVis,chiVis)=shapeGen_VIMS(vertices,vertIndices,
                                                        L_o_n,l_o_n,A_o_n,
                                                        L_s,l_s,A_s,
                                                        dn,
                                                        num.array([oxv,oyv]), ov,ova,
                                                        ps_vis,
                                                        inot,jnot,
                                                        dXV,dYV,dZV,
                                                        lons,lats,
                                                        pTimesV,
                                                        imData[0],vis=True,mask=imData[2],az_adjust=az_adjust)
    (imageIR,poly3d,colours,rot_vert,vertsInPixelIR,chiIR) =shapeGen_VIMS(vertices,vertIndices,
                                                        L_o_n,l_o_n,A_o_n,
                                                        L_s,l_s,A_s,
                                                        dn,
                                                        num.array([oxi,oyi]), ov,ova,
                                                        ps_ir,
                                                        inot,jnot,
                                                        dXI,dYI,dZI,
                                                        lons,lats,
                                                        pTimesI,
                                                        imData[1],az_adjust=az_adjust)
  
    print '  ',(chiVis+chiIR)*0.5,'\n'
    return 0.5*(chiVis+chiIR)


def callShapeGenNoVel(r,vertices,vertIndices,ov,ova,ps_vis,ps_ir,inot,jnot,dXV,dYV,dZV,dXI,dYI,dZI,Ls,ls,pTimesV,pTimesI,imData,az_adjust=-90):
    global steps
    steps+=1
    print steps,

    (L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi)=r

    if abs(oxv)>1200. or abs(oyv)>1200. or oxv==num.inf or oxv==-num.inf or oyv==num.inf or oyv==-num.inf:
        return -num.inf
    if abs(oxi)>1200. or abs(oyi)>1200. or oxi==num.inf or oxi==-num.inf or oyi==num.inf or oyi==-num.inf:
        return -num.inf
    if  dn<0 or ov<0 or ov>5 or ova<0 or ova>360. or dn<0:
        return -num.inf
    if L_o_n>360. or L_o_n<0. or l_o_n>90. or l_o_n<-90. or A_o_n>360. or A_o_n<0. or l_s>90. or l_s<-90. or L_s>360. or L_s<0. or A_s<0 or A_s>360.:
        return -num.inf

    print L_o_n,l_o_n,A_o_n
    print '    ',L_s,l_s,A_s
    print '    ',dn
    print '    ',oxv,oyv,oxi,oyi
    print '    ',ov,ova
    (imageVis,poly3d,colours,rot_vert,vertsInPixelVis,chiVis)=shapeGen_VIMS(vertices,vertIndices,
                                                        L_o_n,l_o_n,A_o_n,
                                                        L_s,l_s,A_s,
                                                        dn,
                                                        num.array([oxv,oyv]), ov,ova,
                                                        ps_vis,
                                                        inot,jnot,
                                                        dXV,dYV,dZV,
                                                        lons,lats,
                                                        pTimesV,
                                                        imData[0],vis=True,mask=imData[2],az_adjust=az_adjust)
    (imageIR,poly3d,colours,rot_vert,vertsInPixelIR,chiIR) =shapeGen_VIMS(vertices,vertIndices,
                                                        L_o_n,l_o_n,A_o_n,
                                                        L_s,l_s,A_s,
                                                        dn,
                                                        num.array([oxi,oyi]), ov,ova,
                                                        ps_ir,
                                                        inot,jnot,
                                                        dXI,dYI,dZI,
                                                        lons,lats,
                                                        pTimesI,
                                                        imData[1],az_adjust=az_adjust)
  
    print '  ',(chiVis+chiIR)*0.5,'\n'
    return 0.5*(chiVis+chiIR)





####data input

imageName='2004163T121836_2004163T192848/cv1465667594_1'
#imageName='2004163T193015_2004164T051726/cv1465680977_5'

#az_hacks
azhacks={}
azhacks['2004163T121836_2004163T192848/cv1465649746_1']=[0,0]

azhacks['2004163T121836_2004163T192848/cv1465669944_1']=[-90,-90]
azhacks['2004163T121836_2004163T192848/cv1465670212_1']=[-90,-30]
azhacks['2004163T121836_2004163T192848/cv1465670650_1']=[-80,-20]
azhacks['2004163T121836_2004163T192848/cv1465671285_1']=[0,-0]
azhacks['2004163T121836_2004163T192848/cv1465671448_1']=[-70,-0]
azhacks['2004163T121836_2004163T192848/cv1465671822_1']=[-60,-60]
azhacks['2004163T121836_2004163T192848/cv1465672161_1']=[-30,-60]
azhacks['2004163T121836_2004163T192848/cv1465673600_1']=[-60,-60]
azhacks['2004163T193015_2004164T051726/cv1465677443_1']=[0,+50]
azhacks['2004163T193015_2004164T051726/cv1465677670_1']=[0,+50]
azhacks['2004163T193015_2004164T051726/cv1465678419_1']=[0,+30]
azhacks['2004163T193015_2004164T051726/cv1465678911_1']=[0,190]
azhacks['2004163T193015_2004164T051726/cv1465679413_1']=[0,+40]
azhacks['2004163T193015_2004164T051726/cv1465679675_1']=[0,+60]
azhacks['2004163T193015_2004164T051726/cv1465679932_1']=[0,+120]
azhacks['2004163T193015_2004164T051726/cv1465680977_2']=[0,+80]
azhacks['2004163T193015_2004164T051726/cv1465680977_5']=[0,+80]
azhacks['2004163T193015_2004164T051726/cv1465678419_1']=[30,70]
try:az_hack_o,az_hack_s=azhacks[imageName]
except: az_hack_o,az_hack_s=0.,0.


if len(sys.argv)>1:
    imageName=sys.argv[1]

print
print 'Working with image:'
print imageName
print


(visDataMin,IRDataMin,oxv,oyv,oxi,oyi,distancenot)=getPars('/data/VIMS/covims_0004/procdata/%s.parFile'%(imageName))
offsetsVis=num.array([oxv,oyv])
offsetsIR=num.array([oxi,oyi])


#code operations
showSurface=True
showSources=True
showNormals=True
plotUpdate=False
alpha=0.7
loadBestPoint=True
doFitsVel=False
doFitsNoVel=False
doRandomVel=False
produceMaps=False
findRoughOffsets=False
exitAfterFit=False
singleFrame=True


if len(sys.argv)>1:
    doRandomVel=False
    loadBestPoint=True
    doFitsVel=False
    exitAfterFit=True
    
with pyf.open('/data/VIMS/covims_0004/procdata/%s_vis_mean_strpd.fits'%(imageName)) as han:
    imDataVis=han[0].data

with pyf.open('/data/VIMS/covims_0004/procdata/%s_vis_mean_strpd_bpmap.fits'%(imageName)) as han:
    maskVis=han[0].data
imDataVis*=maskVis

W=num.where(imDataVis<-10000)
imDataVis[W]=num.median(imDataVis)
maxi=num.max(imDataVis)
W=num.where(imDataVis<maxi)
maxi=num.max(imDataVis[W])
imDataVis=num.clip(imDataVis,num.min(imDataVis),maxi)


w=num.where(imDataVis>visDataMin)
imDataVis[w]=256.
W=num.where(imDataVis<256.)
imDataVis[W]=0.0

with pyf.open('/data/VIMS/covims_0004/procdata/%s_ir_mean.fits'%(imageName)) as han:
    imDataIR=han[0].data


W=num.where(imDataIR<-10000)
imDataIR[W]=num.median(imDataIR)
maxi=num.max(imDataIR)
W=num.where(imDataIR<maxi)
maxi=num.max(imDataIR[W])
imDataIR=num.clip(imDataIR,num.min(imDataIR),maxi)


w=num.where(imDataIR>IRDataMin)
imDataIR[w]=256.
W=num.where(imDataIR<256.)
imDataIR[W]=0.0


if imageName=='2004163T121836_2004163T192848/cv1465651336_1':
    imDataIR[:,0]=0.0


###hack to handle the post-flyby data
if '2004163T193015_2004164T051726' in imageName:
    print 'swap in y'
    imDataVis=imDataVis[::-1,:]
    imDataIR=imDataIR[::-1,:]
    maskVis=maskVis[::-1,:]
elif '2004163T121836_2004163T192848' in imageName:
    print 'swap in x'
    imDataVis=imDataVis[:,::-1]
    imDataIR=imDataIR[:,::-1]
    maskVis=maskVis[:,::-1]

###hack to handle the shadow on Jason
if imageName=='2004163T121836_2004163T192848/cv1465665563_1':
    imDataVis[1:13,10:15]=256.
    imDataIR[9:15,8:12]=256.
imData=num.concatenate([[imDataVis],[imDataIR],[maskVis]])
(dump,A,B)=imData.shape


with open('/data/VIMS/covims_0004/procdata/%s_vis.campt.pickle'%(imageName)) as han:
    latLongObjVis=pick.load(han)

with open('/data/VIMS/covims_0004/procdata/%s_ir.campt.pickle'%(imageName)) as han:
    latLongObjIR=pick.load(han)


#load the shape model

###use the following as reference for plotting
#http://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
if not doFitsVel and not doFitsNoVel:
    x=1
else:
    x=1
n=(64*x)+1
with open('/data/PhoebeShapeModel/CO_SA_ISSNA_5_PHOEBESHAPE_V2_0/data/phoebe_ver%sq.tab'%(n-1)) as han:
    data=han.readlines()

vertices=[]
for i in range(1,6*n**2+1):
    s=data[i].split()
    v=[float(s[1]),float(s[2]),float(s[3])]
    vertices.append(v)

vertices=num.array(vertices)


vertIndices=[]
for i in range(6*n**2+2,len(data)):
    s=data[i].split()
    vertIndices.append([int(float(s[1]))-1,int(float(s[2]))-1,int(float(s[3]))-1])
vertIndices=num.array(vertIndices)

(mids,normals)=midsNormals(vertices[vertIndices])
lonsNotMod=(num.arctan2(mids[:,1],mids[:,0])*r2d)
lons=lonsNotMod%360
lats=num.arcsin(mids[:,2]/(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5)*r2d



sampResolutionsVis=[]
#line is i, sample is j
CassVis=[]
for i in range(len(latLongObjVis)):
    for j in range(len(latLongObjVis[i])):
        if latLongObjVis[i][j]['SpacecraftPosition']<>None:
            sampResolutionsVis.append(latLongObjVis[i][j]['SampleResolution'])
            [X,Y,Z]=latLongObjVis[i][j]['SpacecraftPosition']
            CassVis.append([X,Y,Z,
                            pll.MJD(latLongObjVis[i][j]['UTC']),
                            latLongObjVis[i][j]['TargetCenterDistance'],
                            latLongObjVis[i][j]['Line ']-1,latLongObjVis[i][j]['Sample ']-1])
CassVis=num.array(CassVis)
args=num.argsort(CassVis[:,3])
CassVis=CassVis[args]

sampResolutionsIR=[]
CassIR=[]
for i in range(len(latLongObjIR)):
    for j in range(len(latLongObjIR[i])):
        if latLongObjIR[i][j]['SpacecraftPosition']<>None:
            sampResolutionsIR.append(latLongObjIR[i][j]['SampleResolution'])
            [X,Y,Z]=latLongObjIR[i][j]['SpacecraftPosition']
            CassIR.append([X,Y,Z,
                            pll.MJD(latLongObjIR[i][j]['UTC']),
                            latLongObjIR[i][j]['TargetCenterDistance'],
                            latLongObjIR[i][j]['Line ']-1,latLongObjIR[i][j]['Sample ']-1])

CassIR=num.array(CassIR)
args=num.argsort(CassIR[:,3])
CassIR=CassIR[args]



sampResolutionsIR=num.array(sampResolutionsIR)
sampResolutionsVis=num.array(sampResolutionsVis)
argI=num.argmin(sampResolutionsIR)
timeDelt=num.abs(CassVis[:,3]-CassIR[argI][3])
W=num.where(timeDelt==num.min(timeDelt))

print ' IR nominal sample resolution:',sampResolutionsIR[argI]
print 'Vis nominal sample resolution:',num.min(sampResolutionsVis[W])
resRat=num.min(sampResolutionsVis[W])/sampResolutionsIR[argI]
print resRat


###inot=int(CassIR[0,5])
###jnot=int(CassIR[0,6])
###Tnot=CassIR[0,3]
w=[[0]]
#mid_i=int(num.median(CassIR[:,5]))
#w=num.where(CassIR[:,5]==mid_i)
#mid_j=int(num.median(CassIR[w][:,6]))
#w=num.where((CassIR[:,5]==mid_i)&(CassIR[:,6]==mid_j))
inot=int(CassIR[w][:,5])#same as mid_i
jnot=int(CassIR[w][:,6])#same as mid_j
Tnot=CassIR[w][0][3]



pixelTimesVis=getPixTimes(CassVis,imDataVis.shape[0],imDataVis.shape[1],Tnot,inot,jnot)
pixelTimesIR=getPixTimes(CassIR,imDataIR.shape[0],imDataIR.shape[1],Tnot,inot,jnot)
#for i in range(len(CassIR)):
#    print CassIR[i][3],pixelTimesIR[CassIR[i][5],CassIR[i][6]],CassIR[i][5],CassIR[i][6]
#sys.exit()


#distance units are in km, resolution units in km/pix
print 'Cassini xyz coordinates'
print CassIR[w][0][:3]
[Xnot,Ynot,Znot]=CassIR[w][0][:3]#latLongObjIR[inot][jnot]['SpacecraftPosition']/2
if distancenot==None:
    distancenot=CassIR[w][0][4]#latLongObjIR[inot][jnot]['TargetCenterDistance']/2
else:
    Xnot*=distancenot/CassIR[w][0][0]
    Ynot*=distancenot/CassIR[w][0][0]
    Znot*=distancenot/CassIR[w][0][0]
#sampleResolutionnot=latLongObjIR[inot][jnot]['SampleResolution']/1000.

#nominal values from the VIMS manual are 0.17*3 and 0.5 mrad
if resRat<0.5:
    pixScaleVis=0.17#*.33/.35
else:
    pixScaleVis=0.17*3
pixScaleIR=0.5
#sampleResolutionIRnot=irPixScale*distancenot/1000.
#sampleResolutionVisnot=visPixScale*distancenot/1000.

(junk,deltaXVis,deltaYVis,deltaZVis)=getVels(CassVis,pixelTimesVis,Xnot,Ynot,Znot,inot,jnot,polyOrder=2)
(spaceCraftVectornot,deltaXIR,deltaYIR,deltaZIR)=getVels(CassIR,pixelTimesIR,Xnot,Ynot,Znot,inot,jnot,polyOrder=2)


long_o_not=latLongObjIR[inot][jnot]['SubSpacecraftLongitude']
lat_o_not=latLongObjIR[inot][jnot]['SubSpacecraftLatitude']
az_o_not=latLongObjIR[inot][jnot]['SpacecraftAzimuth']+az_hack_o
long_s=latLongObjIR[inot][jnot]['SubSolarLongitude']
lat_s=latLongObjIR[inot][jnot]['SubSolarLatitude']
az_s=latLongObjIR[inot][jnot]['SubSolarAzimuth']+az_hack_s


#velocities
#offsetsVis=num.array([104.,291.]) #shift in vertical, shift in horizontal
#offsetsIR=num.array([112.,167.]) #shift in vertical, shift in horizontal
offsetVel=0.0 #km/s
#0.0 is an expansion in verticle, 180 is a contraction in vertical
#270 is an expansion in horizontal, 90 is a contraction in horizontal
offsetVelAngle=270.0 #avoid starting this near 0.

#hacks
if az_s==0.:
    az_s=168.
if az_o_not==0:
    az_o_not=44.
az_s=az_s%360.
az_o_not=az_o_not%360

print "Spice Kernel based angles:"
print long_o_not,lat_o_not,az_o_not
print long_s,lat_s,az_s
print

#sys.exit()



if findRoughOffsets:
    best=[-1000000.,-1,-1]
    offXI=num.linspace(-700,700,21)
    offYI=num.linspace(-700,700,21)
    for i in range(len(offXI)):
        for j in range(len(offYI)):
            (imageIR,poly3d,colours,rot_vert,vertsInPixelIR,chi)=shapeGen_VIMS(vertices,vertIndices,
                                                    long_o_not,lat_o_not,az_o_not,
                                                    long_s,lat_s,az_s,
                                                    distancenot,
                                                    num.array([offXI[i],offYI[j]]), offsetVel,offsetVelAngle,
                                                    pixScaleIR,
                                                    inot,jnot,
                                                    deltaXIR,deltaYIR,deltaZIR,
                                                    lons,lats,
                                                    pixelTimesVis,
                                                    imData[0],vis=False)
            print offXI[i],offYI[j],chi
            if chi>best[0]:
                best=[chi,offXI[i],offYI[j]]
    sys.exit()

#free pars are:
#spacecraft long,lat,az
#solar long, lat, az
#spacecraft distancenot
#offsets (x,y)
#offset velocity
#offset angle
#total 11 free pars

if doFitsVel:
    x,y=0.,0.
    
    nDim=13
    nWalkers=nDim*4
    nBurn=100
    nStep=100

    steps=0
    offWidth=20.
    angWidth=10.
    velWidth=0.1
    distanceWidth=distancenot*0.1
    
    r0=[]
    widthArr=num.array([angWidth,angWidth,angWidth,angWidth,angWidth,angWidth,
                        distanceWidth,offWidth,offWidth,offWidth,offWidth,velWidth,angWidth*10]) #when eyeballing at ~270 degrees, best fit was -38
    for i in range(nWalkers):
        entry=num.array([long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,
                         distancenot,offsetsVis[0],offsetsVis[1],offsetsIR[0],offsetsIR[1],offsetVel,offsetVelAngle])
        entry+=sci.randn(nDim)*widthArr
        entry[len(entry)-2]=abs(entry[len(entry)-2])
        entry[2]=entry[2]%360
        entry[5]=entry[5]%360
        entry[len(entry)-1]=entry[len(entry)-1]%360. #making sure the offset velocity angle is between 0 and 360
        if entry[len(entry)-1]>360:
            sys.exit()
        r0.append(entry)
    r0=num.array(r0)

        
    sampler = emcee.EnsembleSampler(nWalkers, nDim, callShapeGen, args=[vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData,-90])

    pos, prob, state = sampler.run_mcmc(r0, nBurn)

    sampler.reset()

    pos, prob, state = sampler.run_mcmc(pos, nStep, rstate0=state)

    samps=sampler.chain
    probs=sampler.lnprobability

    with open('/data/VIMS/covims_0004/procdata/%s.fit_pickle'%(imageName),'w+') as outHan:
        pick.dump([samps,probs],outHan)
    
    print sampler.chain
    print
    print sampler.lnprobability

    (bestPoint,goodSamps)=getFit('/data/VIMS/covims_0004/procdata/%s.fit_pickle'%(imageName))
    print bestPoint
    if exitAfterFit: sys.exit()

if doFitsNoVel:
    x,y=0.,0.
    
    nDim=11
    nWalkers=nDim*4
    nBurn=30
    nStep=30

    steps=0
    offWidth=20.
    angWidth=10.
    distanceWidth=distancenot*0.1
    
    r0=[]
    widthArr=num.array([angWidth,angWidth,angWidth,angWidth,angWidth,angWidth,
                        distanceWidth,offWidth,offWidth,offWidth,offWidth])
    for i in range(nWalkers):
        entry=num.array([long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,
                         distancenot,offsetsVis[0],offsetsVis[1],offsetsIR[0],offsetsIR[1]])
        entry+=sci.randn(nDim)*widthArr
        #entry[10]=entry[10]%360. #making sure the offset velocity angle is between 0 and 360
        r0.append(entry)
    r0=num.array(r0)
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, callShapeGenNoVel, args=[vertices,vertIndices,offsetVel,offsetVelAngle,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData,-90])

    pos, prob, state = sampler.run_mcmc(r0, nBurn)

    sampler.reset()

    pos, prob, state = sampler.run_mcmc(pos, nStep, rstate0=state)

    samps=sampler.chain
    probs=sampler.lnprobability

    with open('/data/VIMS/covims_0004/procdata/%s.noVel_fit_pickle'%(imageName),'w+') as outHan:
        pick.dump([samps,probs],outHan)
    
    print sampler.chain
    print
    print sampler.lnprobability

    (bestPoint,goodSamps)=getFit('/data/VIMS/covims_0004/procdata/%s.noVel_fit_pickle'%(imageName))
    print bestPoint






###modelling free parameters are:
# space craft distance
# subspacecraft long,lat,az
# subsolar long, lat, az
# y,z image offsets in km
# offsetVel in km/s
# offsetVelAngle in degrees
###Chosen parameters are:
# polyOrder for the velocity modelling


if doRandomVel:
    print 'Doing random sampling'
    nRand=8000
    steps=0
    (bestPoint,goodSamps)=getFit('/data/VIMS/covims_0004/procdata/%s.fit_pickle'%(imageName))
    bestPoint[len(bestPoint)-1]=callShapeGen(bestPoint[:len(bestPoint)-1],vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData,-90)
    (ll,l)=goodSamps.shape
    print 'Finding better fits than '+str(bestPoint[len(bestPoint)-1])
    steps=0
    
    r0=[]
    widthArr=bestPoint[:len(bestPoint)-1]-([long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offsetsVis[0],offsetsVis[1],offsetsIR[0],offsetsIR[1],offsetVel,offsetVelAngle])

    for i in range(nRand):
        entry=bestPoint[:len(bestPoint)-1]+sci.randn(len(widthArr))*widthArr/4.
        entry[len(entry)-1]=entry[len(entry)-1]%360. #making sure the offset velocity angle is between 0 and 360
        r0.append(entry)
    r0=num.array(r0)


    for i in range(nRand):
        chi=callShapeGen(r0[i],vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData,-90)
        if chi>=bestPoint[len(bestPoint)-1]:
            with open('/data/VIMS/covims_0004/procdata/%s.rand'%(imageName),'a+') as handle:
                print >>handle,r0[i],chi
                print '****'
                print
    sys.exit()

if loadBestPoint:# and imageName<>'2004163T193015_2004164T051726/cv1465680977_5':

    steps=0
    (bestPoint,goodSamps)=getFit('/data/VIMS/covims_0004/procdata/%s.fit_pickle'%(imageName))
    print bestPoint
    [long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offXV,offYV,offXI,offYI,offsetVel,offsetVelAngle,chi]=bestPoint
    chi=callShapeGen(bestPoint[:len(bestPoint)-1],vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData,-90)
    print 'EMCEE best chi:',chi

    if os.path.isfile('/data/VIMS/covims_0004/procdata/%s.rand'%(imageName)):
        with open('/data/VIMS/covims_0004/procdata/%s.rand'%(imageName)) as randHan:
            randoms=randHan.readlines()
        for i in range(0,len(randoms),4):
            s=''
            for j in range(4):
                s+=randoms[i+j].split('\n')[0]
            s=s.split()
            rchi=float(s[len(s)-1])
            #print rchi
            if rchi>chi:
                x=[]
                for j in range(1,len(s)-1):
                    x.append(float(s[j].split(']')[0]))
                x=num.array(x)
                [long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offXV,offYV,offXI,offYI,offsetVel,offsetVelAngle]=x
                chi=rchi
        print 'RANDOM best chi:',chi
        print
    offsetsVis=num.array([offXV,offYV])
    offsetsIR =num.array([offXI,offYI])
    
    



(imageVis,poly3d_r,colours,rot_vert_r,crap,chi)=shapeGen_VIMS(vertices,vertIndices,
                                                    long_o_not,lat_o_not,az_o_not,
                                                    long_s,lat_s,az_s,
                                                    distancenot,
                                                    offsetsVis, offsetVel,offsetVelAngle,
                                                    pixScaleVis,
                                                    inot,jnot,
                                                    deltaXVis,deltaYVis,deltaZVis,
                                                    lons,lats,
                                                    pixelTimesVis,
                                                    imData[0],vis=True,mask=imData[2],az_adjust=-0)
(imageVis,poly3d,colours,rot_vert,vertsInPixelVis,chi)=shapeGen_VIMS(vertices,vertIndices,
                                                    long_o_not,lat_o_not,az_o_not,
                                                    long_s,lat_s,az_s,
                                                    distancenot,
                                                    offsetsVis, offsetVel,offsetVelAngle,
                                                    pixScaleVis,
                                                    inot,jnot,
                                                    deltaXVis,deltaYVis,deltaZVis,
                                                    lons,lats,
                                                    pixelTimesVis,
                                                    imData[0],vis=True,mask=imData[2],az_adjust=-90)
(imageIR,poly3d,colours,rot_vert,vertsInPixelIR,chi)=shapeGen_VIMS(vertices,vertIndices,
                                                    long_o_not,lat_o_not,az_o_not,
                                                    long_s,lat_s,az_s,
                                                    distancenot,
                                                    offsetsIR, offsetVel,offsetVelAngle,
                                                    pixScaleIR,
                                                    inot,jnot,
                                                    deltaXIR,deltaYIR,deltaZIR,
                                                    lons,lats,
                                                    pixelTimesVis,
                                                    imData[0],vis=False,az_adjust=-90)
steps=0
callShapeGen((long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offsetsVis[0],offsetsVis[1],offsetsIR[0],offsetsIR[1],offsetVel,offsetVelAngle),vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData,-90)




if produceMaps :#and loadBestPoint:
    with pyf.open('/data/VIMS/covims_0004/procdata/'+imageName+'_ir.fits') as shan:
        specData=shan[0].data
    if '2004163T193015_2004164T051726' in imageName:
        print 'not swaitch'
        #specData=specData[:,:,:]
    elif '2004163T121836_2004163T192848' in imageName:
        print 'swaitch'
        specData=specData[:,::-1,::-1]


    (channels,A,B)=specData.shape

    waterDepth=[]
    for l in range(A):
        waterDepth.append([])
        for s in range(B):
            (w,spec)=specAnalysis.getSpec(specData,l=l,s=s)
            med=num.nanmedian(spec)
            #print l,s,med,imData[1,l,s]

            if med>0.0005:
                (h,junk)=specAnalysis.water(w,spec)
                if h>0.0:
                    waterDepth[len(waterDepth)-1].append(h)
                else: waterDepth[len(waterDepth)-1].append(0.0)
            else:
                waterDepth[len(waterDepth)-1].append(0.0)
    waterDepth=num.array(waterDepth)

    """
    for i in range(len(waterDepth)):
        for j in range(len(waterDepth[i])):
            if len(vertsInPixelIR[i][j])>0:
                print i,j
    print
    for i in range(len(waterDepth)):
        for j in range(len(waterDepth[i])):
            if waterDepth[i,j]<>0:
                print i,j
    #sys.exit()
    """

    x=[]
    y=[]
    for i in range(len(vertsInPixelIR)):
        for j in range(len(vertsInPixelIR[i])):
            if len(vertsInPixelIR[i][j])>0 and waterDepth[i,j]>0:

                m=mids[vertsInPixelIR[i][j].astype('int')]
                h=(m[:,0]**2+m[:,1]**2+m[:,2]**2)**0.5
                hmed=num.median(h)
                lon=num.median(lons[vertsInPixelIR[i][j].astype('int')])
                lat=num.median(lats[vertsInPixelIR[i][j].astype('int')])

                #if not num.isnan(hmed): print len(h),len(vertsInPixelIR[i][j])
                #    if not num.isnan(hmed) and waterDepth[i,j]<>0:
                x.append(hmed)
                y.append(waterDepth[i,j])

                for k in range(len(vertsInPixelIR[i][j])):
                    #element0==1 we get bright red
                    #element1==1 we get bright green
                    #element2==1 we get bright navy blue
                    if waterDepth[i,j]<0.22:
                        colours[vertsInPixelIR[i][j][k]][0]=1.0 #red
                    elif waterDepth[i,j]<0.32:
                        colours[vertsInPixelIR[i][j][k]][1]=1.0 #green
                    else:
                        colours[vertsInPixelIR[i][j][k]][2]=1.0 #blue


    fig2=pyl.figure(2)    
    pyl.scatter(x,y)
    pyl.xlabel('Depth (km)')
    pyl.ylabel('Water Absorption Depth')



if singleFrame:
    #single window plot
    fig=pyl.figure(1,figsize=(10,10))
    pyl.title(imageName+'\n Purple=Image Data, Orange=Model')
    fig.subplots_adjust(wspace=0,hspace=0)
    ax1=fig.add_subplot(1,1,1,projection='3d')
    #if showSurface:
    #    ax1.set_axis_bgcolor('0.0')
    #    ax1.set_axis_off()
    ax1.set_aspect('equal')
    ax1.set_xlim(-120,120)
    ax1.set_ylim(-120,120)
    ax1.set_zlim(-120,120)
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    azim=0.0
    elev=0.0
    ax1.view_init(azim=azim, elev=elev)
    collection=Poly3DCollection(poly3d_r,linewidths=0.0,facecolors=colours)
    collection.set_alpha(1.0)
    pylCollection=ax1.add_collection3d(collection)

    pyl.show()
    sys.exit()


#4 panel plot!
fig=pyl.figure(1,figsize=(10,10))
pyl.title(imageName+'\n Purple=Image Data, Orange=Model')
fig.subplots_adjust(wspace=0,hspace=0)
ax1=fig.add_subplot(2,2,1,projection='3d')
ax2=fig.add_subplot(2,2,2,projection='3d')
ax3=fig.add_subplot(2,2,3)
ax4=fig.add_subplot(2,2,4)
#fig.canvas.mpl_connect('key_press_event',event_handler)


if showSurface:
    ax1.set_axis_bgcolor('0.0')
    ax1.set_axis_off()
ax1.set_aspect('equal')
ax2.set_aspect('equal')
ax3.set_aspect('equal')
ax4.set_aspect('equal')

ax1.set_xlim(-120,120)
ax1.set_ylim(-120,120)
ax1.set_zlim(-120,120)
ax2.set_xlim(-120,120)
ax2.set_ylim(-120,120)
ax2.set_zlim(-120,120)

azim=0.0
elev=0.0
ax1.view_init(azim=azim, elev=elev)
ax2.view_init(azim=azim, elev=elev)

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')



collection=Poly3DCollection(poly3d_r,linewidths=0.0,facecolors=colours)
collection.set_alpha(alpha)
pylCollection=ax1.add_collection3d(collection)

if showSurface:
    vertSpacing=2#0
else:
    vertSpacing=1
ax2.scatter(rot_vert_r[:,0][::vertSpacing],
            rot_vert_r[:,1][::vertSpacing],
            rot_vert_r[:,2][::vertSpacing],
            marker='.',s=2,c='k',alpha=0.5)



plotData=imData[0]*0.0
w=num.where((imData[0]>0)&(imageVis>0))
plotData[w]=1.0
if len(w[0])==0:
    plotData[0,0]=1.0
    plotData[A-1,0]=1.0
    plotData[0,B-1]=1.0
    plotData[A-1,B-1]=1.0
w=num.where((imData[0]>0)&(imageVis==0))
plotData[w]=0.3
w=num.where((imData[0]==0)&(imageVis>0))
plotData[w]=0.7
realImPlot=ax3.imshow(plotData,aspect='equal',interpolation='nearest',origin='lower')
realImPlot.set_cmap('CMRmap')
ax3.set_xlabel('Visual')

plotData=imData[1]*0.0
w=num.where((imData[1]>0)&(imageIR>0))
plotData[w]=1.0
if len(w[0])==0:
    plotData[0,0]=1.0
    plotData[A-1,0]=1.0
    plotData[0,B-1]=1.0
    plotData[A-1,B-1]=1.0
w=num.where((imData[1]>0)&(imageIR==0))
plotData[w]=0.3
w=num.where((imData[1]==0)&(imageIR>0))
plotData[w]=0.7
modelImPlot=ax4.imshow(plotData,aspect='equal',interpolation='nearest',origin='lower')
modelImPlot.set_cmap('CMRmap')
ax4.set_xlabel('IR')


pyl.show()
             
