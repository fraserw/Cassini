#! /usr/bin/env python

import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyl
from matplotlib import backend_bases
import numpy as num, scipy as sci
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys
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
import iminuit
#from numpy import linalg
from iminuit import Minuit

class dataObject():
    def __init__(self,vertices,vertIndices,ps_vis,ps_ir,inot,jnot,dXV,dYV,dZV,dXI,dYI,dZI,Ls,ls,pTimesV,pTimesI,imData):
        self.vertices=vertices
        self.vertIndices=vertIndices
        self.ps_vis=ps_vis
        self.ps_ir=ps_ir
        self.inot=inot
        self.jnot=jnot
        self.dXV=dXV
        self.dYV=dYV
        self.dZV=dZV
        self.dXI=dXI
        self.dYI=dYI
        self.dZI=dZI
        self.Ls=Ls,
        self.ls=ls
        self.pTimesV=pTimesV
        self.pTimesI=pTimesI
        self.imData=imData
        self.chi=None
        
    def callShapeGen(self,L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi,ov,ova):
        r=(L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi,ov,ova)
        chi=callShapeGen(r,self.vertices,self.vertIndices,self.ps_vis,self.ps_ir,self.inot,self.jnot,self.dXV,self.dYV,self.dZV,self.dXI,self.dYI,self.dZI,self.Ls,self.ls,self.pTimesV,self.pTimesI,self.imData)
        self.chi=chi
        return -chi

    def callShapeGenNoVel(self,L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi):
        r=(L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi)
        chi=callShapeGenNoVel(r,self.vertices,self.vertIndices,0.0,270.,self.ps_vis,self.ps_ir,self.inot,self.jnot,self.dXV,self.dYV,self.dZV,self.dXI,self.dYI,self.dZI,self.Ls,self.ls,self.pTimesV,self.pTimesI,self.imData)
        self.chi=chi
        return -chi
    
def callShapeGen(r,vertices,vertIndices,ps_vis,ps_ir,inot,jnot,dXV,dYV,dZV,dXI,dYI,dZI,Ls,ls,pTimesV,pTimesI,imData):
    global steps
    steps+=1
    print steps,

    (L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi,ov,ova)=r

    if abs(oxv)>700. or abs(oyv)>700. or oxv==num.inf or oxv==-num.inf or oyv==num.inf or oyv==-num.inf:
        return -num.inf
    if abs(oxi)>700. or abs(oyi)>700. or oxi==num.inf or oxi==-num.inf or oyi==num.inf or oyi==-num.inf:
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
                                                        imData[0],vis=True)
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
                                                        imData[1])
  
    print '  ',(chiVis+chiIR)*0.5,'\n'
    return 0.5*(chiVis+chiIR)


def callShapeGenNoVel(r,vertices,vertIndices,ov,ova,ps_vis,ps_ir,inot,jnot,dXV,dYV,dZV,dXI,dYI,dZI,Ls,ls,pTimesV,pTimesI,imData):
    global steps
    steps+=1
    print steps,

    (L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,oxv,oyv,oxi,oyi)=r

    if abs(oxv)>700. or abs(oyv)>700. or oxv==num.inf or oxv==-num.inf or oyv==num.inf or oyv==-num.inf:
        return -num.inf
    if abs(oxi)>700. or abs(oyi)>700. or oxi==num.inf or oxi==-num.inf or oyi==num.inf or oyi==-num.inf:
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
                                                        imData[0],vis=True)
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
                                                        imData[1])
  
    print '  ',(chiVis+chiIR)*0.5,'\n'
    return 0.5*(chiVis+chiIR)





####data input parameters
#imageName='2004163T121836_2004163T192848/cv1465671822_1'
#visDataMin=0.001
#IRDataMin=0.001
#weird orientation for the below image
#imageName='2004163T121836_2004163T192848/cv1465669068_1'
#visDataMin=0.022
#IRDataMin=0.002
imageName='2004163T121836_2004163T192848/cv1465670212_1'
visDataMin=0.004
IRDataMin=0.004

#code operations
showSurface=True
showSources=True
showNormals=True
plotUpdate=False
alpha=0.2
loadBestPoint=True
doFitsVel=False
doFitsNoVel=False




with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_vis_mean.fits'%(imageName)) as han:
    imDataVis=han[0].data

with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_vis_mean_bpmap.fits'%(imageName)) as han:
    mask=han[0].data
imDataVis*=mask

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

with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_ir_mean.fits'%(imageName)) as han:
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

imData=num.concatenate([[imDataVis],[imDataIR]])
(dump,A,B)=imData.shape


with open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_vis.campt.pickle'%(imageName)) as han:
    latLongObjVis=pick.load(han)

with open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_ir.campt.pickle'%(imageName)) as han:
    latLongObjIR=pick.load(han)


#load the shape model

###use the following as reference for plotting
#http://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
if not doFitsVel and not doFitsNoVel:
    x=1
else:
    x=1
n=(64*x)+1
with open('/Volumes/data/PhoebeShapeModel/CO_SA_ISSNA_5_PHOEBESHAPE_V2_0/data/phoebe_ver%sq.tab'%(n-1)) as han:
    data=han.readlines()

vertices=[]
for i in range(1,6*n**2+1):
    s=data[i].split()
    v=[float(s[1]),float(s[2]),float(s[3])]
    vertices.append(v)

vertices=num.array(vertices)


#get the vertices indices
vertIndices=[]
for i in range(6*n**2+2,len(data)):
    s=data[i].split()
    vertIndices.append([int(float(s[1]))-1,int(float(s[2]))-1,int(float(s[3]))-1])
vertIndices=num.array(vertIndices)

(mids,normals)=midsNormals(vertices[vertIndices])
lons=(num.arctan2(mids[:,1],mids[:,0])*r2d)%360
lats=num.arcsin(mids[:,2]/(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5)*r2d




#plot!
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
pixelTimesVis=getPixTimes(CassVis,A,B,Tnot,inot,jnot)
pixelTimesIR=getPixTimes(CassIR,A,B,Tnot,inot,jnot)



#distance units are in km, resolution units in km/pix
print CassIR[w][0][:3]
[Xnot,Ynot,Znot]=CassIR[w][0][:3]#latLongObjIR[inot][jnot]['SpacecraftPosition']/2
distancenot=CassIR[w][0][0]#latLongObjIR[inot][jnot]['TargetCenterDistance']/2
#sampleResolutionnot=latLongObjIR[inot][jnot]['SampleResolution']/1000.

#nominal values from the VIMS manual are 0.17*3 and 0.5 mrad
if resRat<0.5:
    pixScaleVis=0.17
else:
    pixScaleVis=0.17*3
pixScaleIR=0.5
#sampleResolutionIRnot=irPixScale*distancenot/1000.
#sampleResolutionVisnot=visPixScale*distancenot/1000.

(junk,deltaXVis,deltaYVis,deltaZVis)=getVels(CassVis,pixelTimesVis,Xnot,Ynot,Znot,inot,jnot,polyOrder=2)
(spaceCraftVectornot,deltaXIR,deltaYIR,deltaZIR)=getVels(CassIR,pixelTimesIR,Xnot,Ynot,Znot,inot,jnot,polyOrder=2)


long_o_not=latLongObjIR[inot][jnot]['SubSpacecraftLongitude']
lat_o_not=latLongObjIR[inot][jnot]['SubSpacecraftLatitude']
az_o_not=latLongObjIR[inot][jnot]['SpacecraftAzimuth']#+120.
long_s=latLongObjIR[inot][jnot]['SubSolarLongitude']
lat_s=latLongObjIR[inot][jnot]['SubSolarLatitude']
az_s=latLongObjIR[inot][jnot]['SubSolarAzimuth']#+120.

offsetsVis=num.array([104.,291.]) #shift in vertical, shift in horizontal
offsetsIR=num.array([112.,167.]) #shift in vertical, shift in horizontal
offsetVel=0.0 #km/s
#0.0 is an expansion in verticle, 180 is a contraction in vertical
#270 is an expansion in horizontal, 90 is a contraction in horizontal
offsetVelAngle=270.0 #avoid starting this near 0.


#free pars are:
#spacecraft long,lat,az
#solar long, lat, az
#spacecraft distancenot
#offsets (x,y)
#offset velocity
#offset angle
#total 11 free pars
    





###modelling free parameters are:
# space craft distance
# subspacecraft long,lat,az
# subsolar long, lat, az
# y,z image offsets in km
# offsetVel in km/s
# offsetVelAngle in degrees
###Chosen parameters are:
# polyOrder for the velocity modelling

if loadBestPoint:
    (bestPoint,goodSamps)=getFit('/Users/fraserw/data/VIMS/covims_0004/procdata/%s.fit_pickle'%(imageName))
    (ll,l)=goodSamps.shape
    #w=num.where(goodSamps[:,13]>bestPoint[13]-2)
    #for i in range(len(w[0])):
    #    print goodSamps[w[0][i]]
    #sys.exit()
    [long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offXV,offYV,offXI,offYI,offsetVel,offsetVelAngle,chi]=bestPoint
    offsetsVis=num.array([offXV,offYV])
    offsetsIR =num.array([offXI,offYI])
    

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
                                                    imData[0],vis=True)
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
                                                    imData[0],vis=False)
steps=0

#callShapeGen((long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offsetsVis[0],offsetsVis[1],offsetsIR[0],offsetsIR[1],offsetVel,offsetVelAngle),vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData)




minObject=dataObject(vertices,vertIndices,pixScaleVis,pixScaleIR,inot,jnot,deltaXVis,deltaYVis,deltaZVis,deltaXIR,deltaYIR,deltaZIR,lons,lats,pixelTimesVis,pixelTimesIR,imData)
#minObject.callShapeGenNoVel(long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offsetsVis[0],offsetsVis[1],offsetsIR[0],offsetsIR[1])
angleError=1000.
distanceError=5000.
offsetError=5000.
velError=100.
minuit=Minuit(minObject.callShapeGen,
              L_o_n=long_o_not,
              l_o_n=lat_o_not,
              A_o_n=az_o_not,
              L_s=long_s,
              l_s=lat_s,
              A_s=az_s,
              dn=distancenot,
              oxv=offsetsVis[0],
              oyv=offsetsVis[1],
              oxi=offsetsIR[0],
              oyi=offsetsIR[1],
              ov=offsetVel,
              ova=offsetVelAngle,
              error_L_o_n=angleError,
              error_l_o_n=angleError,
              error_A_o_n=angleError,
              error_L_s=angleError,
              error_l_s=angleError,
              error_A_s=angleError,
              error_dn=distanceError,
              error_oxv=offsetError,
              error_oyv=offsetError,
              error_oxi=offsetError,
              error_oyi=offsetError,
              error_ov=velError,
              error_ova=angleError,
              limit_L_o_n=(0.,360.),
              limit_l_o_n=(-90.,90.),
              limit_A_o_n=(0.,360.),
              limit_L_s=(0.,360.),
              limit_l_s=(-90.,90.),
              limit_A_s=(0.,360.),
              limit_ov=(0.,0.4),
              limit_ova=(0.,360.))
            
              
minuit.migrad()
print minObject.chi
sys.exit()







 

collection=Poly3DCollection(poly3d,linewidths=0.0,facecolors=colours)
collection.set_alpha(alpha)
pylCollection=ax1.add_collection3d(collection)

if showSurface:
    vertSpacing=2#0
else:
    vertSpacing=1
ax2.scatter(rot_vert[:,0][::vertSpacing],
            rot_vert[:,1][::vertSpacing],
            rot_vert[:,2][::vertSpacing],
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
realImPlot=ax3.imshow(plotData[::-1,:],aspect='equal',interpolation='nearest')
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
modelImPlot=ax4.imshow(plotData[::-1,:],aspect='equal',interpolation='nearest')
modelImPlot.set_cmap('CMRmap')
ax4.set_xlabel('IR')


pyl.show()
             
