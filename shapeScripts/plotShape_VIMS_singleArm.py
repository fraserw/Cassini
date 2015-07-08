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

from shapeLib import *
from lineFitUtils import poly
#from numpy import linalg

def callShapeGen(r,vertices,vertIndices,ov,ova,sr_n,inot,jnot,dX,dY,dZ,Ls,ls,pTimes,imData):
    global steps

    (L_o_n,l_o_n,A_o_n,L_s,l_s,A_s,dn,ox,oy)=r

    steps+=1
    if abs(ox)>500. or abs(oy)>500. or ox==num.inf or ox==-num.inf or oy==num.inf or oy==-num.inf or dn<0 or ov<0 or ov>5 or ova<0 or ova>360. or sr_n<0:
        return -num.inf
    if L_o_n>360. or L_o_n<0. or l_o_n>90. or l_o_n<-90. or A_o_n>360. or A_o_n<0. or l_s>90. or l_s<-90. or L_s>360. or L_s<0. or A_s<0 or A_s>360.:
        return -num.inf

    print steps,L_o_n,l_o_n,A_o_n
    print '    ',L_s,l_s,A_s
    print '    ',dn,ox,oy
    print '    ',ov,ova,sr_n
    (image,poly3d,colours,rot_vert,chi)=shapeGen_VIMS(vertices,vertIndices,
                                      L_o_n,l_o_n,A_o_n,
                                      L_s,l_s,A_s,
                                      dn,
                                      num.array([ox,oy]), ov,ova,
                                      sr_n,
                                      inot,jnot,
                                      deltaX,deltaY,deltaZ,
                                      lons,lats,
                                      pixelTimes,
                                      imData)
  
    print '  ',chi,'\n'
    return 0.5*chi

            
def event_handler(event):
    global vertics,vertIndices,showSurface,showNormals,showSources
    global long_o,lat_o,a_o,long_s,lat_s,a_s
    global sampleResolution
    global imData
    global plotUpdate
    global alpha
    global azim,elev

    print event.key


    if event.key in ['shift','ctrl+','ctrl+shift','control']: return None
    elif event.key=='q': sys.exit()
        
    if event.key=='z':long_o+=15.0
    elif event.key=='Z':long_o-=15.0
    if event.key=='x':lat_o+=15.0
    elif event.key=='X':lat_o-=15.0
    if event.key=='c':a_o+=15.0
    elif event.key=='C':a_o-=15.0

    if event.key=='ctrl+z':long_s+=5.0
    elif event.key=='ctrl+Z':long_s-=5.0
    if event.key=='ctrl+x':lat_s+=5.0
    elif event.key=='ctrl+X':lat_s-=5.0
    if event.key=='ctrl+c':a_s+=5.0
    elif event.key=='ctrl+C':a_s-=5.0

    print long_o%360,lat_o,a_o%360
    print long_s%360.,lat_s,a_s%360.

    
    if event.key=='r':
        long_s+=50.*(random.random()-0.5)
        lat_s+=50.*(random.random()-0.5)
        a_s+=50.*(random.random()-0.5)
        long_o+=50.*(random.random()-0.5)
        lat_o+=50.*(random.random()-0.5)
        a_o+=50.*(random.random()-0.5)
    
        
    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)=shapeGen_VIMS(vertices,vertIndices,
                                                                                     long_o=long_o,lat_o=lat_o,az_o=a_o,
                                                                                     long_s=long_s,lat_s=lat_s,
                                                                                     imData=imData*1.0)

    collection=Poly3DCollection(p3d,linewidths=0.0,facecolors=c3d)
    collection.set_alpha(alpha)
    
    ax1.cla()
    ax2.cla()
    ax3.cla()
    ax4.cla()
    
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
    

    imgPlot=ax3.imshow(256.-image,aspect='equal',interpolation='nearest')
    imgPlot.set_cmap('Greys')
    oImData=imData*0.7
    oImData+=image*0.3
    imgDataPlot=ax4.imshow(256.-oImData,aspect='equal',interpolation='nearest')
    imgDataPlot.set_cmap('Greys')


    fig.canvas.set_window_title(chi)

    
    #w=num.where(c2d==0)
    #ax3.scatter(m2d[w,1],m2d[w,2],marker='.')
    #w=num.where(c2d>0)
    #ax3.scatter(m2d[w,1],m2d[w,2],marker='.',c='r')
    
    if showSurface:
        pylCollection=ax1.add_collection3d(collection)

    if showSurface:
        vertSpacing=2#0
    else:
        vertSpacing=1
    ax2.scatter(rot_vertices[:,0][::vertSpacing],
                rot_vertices[:,1][::vertSpacing],
                rot_vertices[:,2][::vertSpacing],
                marker='.',s=2,c='k',alpha=0.5)

    if showSources:
        ax2.scatter(rot_n_sun[0]*200.,rot_n_sun[1]*200.,rot_n_sun[2]*200.,marker=(5,0,0),c='y',s=70,zorder=0)
        ax2.scatter(n_obs[0]*200.,n_obs[1]*200.,n_obs[2]*200.,marker=(10,0,0),c='g',s=70,zorder=0)
        #ax2.scatter(0.0,0.0,0.0,marker=(10,0,0),c='r',s=50,zorder=0)
    
    if showNormals:
        for i in range(567,568):#0,len(n),45):
            #if angs_s[i]<=num.pi/2.:
            Y=num.concatenate([m[i],m[i]+n[i]*50.]).reshape(2,3)
            ax2.plot(Y[:,0],Y[:,1],Y[:,2],c='r',lw=3)
    ax1.view_init(azim=azim, elev=elev)
    ax2.view_init(azim=azim, elev=elev)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')

    pyl.draw()
    
image='2004163T121836_2004163T192848/cv1465671822_1_vis'

with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_mean.fits'%(image)) as han:
    imData=han[0].data

if 'vis' in image:
    with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_mean_bpmap.fits'%(image)) as han:
        mask=han[0].data
    imData*=mask

W=num.where(imData<-10000)
imData[W]=num.median(imData)
maxi=num.max(imData)
W=num.where(imData<maxi)
maxi=num.max(imData[W])
imData=num.clip(imData,num.min(imData),maxi)


w=num.where(imData>0.001)
imData[w]=256.
W=num.where(imData<256.)
imData[W]=0.0



with open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s.campt.pickle'%(image)) as han:
    latLongObj=pick.load(han)

showSurface=True
showSources=True
showNormals=True
plotUpdate=False
alpha=0.2
doFits=False

#load the shape model

###use the following as reference for plotting
#http://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
if not doFits:
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
fig.subplots_adjust(wspace=0,hspace=0)
ax1=fig.add_subplot(2,2,1,projection='3d')
ax2=fig.add_subplot(2,2,2,projection='3d')
ax3=fig.add_subplot(2,2,3)
ax4=fig.add_subplot(2,2,4)
fig.canvas.mpl_connect('key_press_event',event_handler)


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


#line is i, sample is j
Cass=[]
for i in range(len(latLongObj)):
    for j in range(len(latLongObj[i])):
        if latLongObj[i][j]['SpacecraftPosition']<>None:
            [X,Y,Z]=latLongObj[i][j]['SpacecraftPosition']
            Cass.append([X,Y,Z,
                         pll.MJD(latLongObj[i][j]['UTC']),
                         latLongObj[i][j]['TargetCenterDistance'],
                         latLongObj[i][j]['Line ']-1,latLongObj[i][j]['Sample ']-1])
Cass=num.array(Cass)
args=num.argsort(Cass[:,3])
Cass=Cass[args]

#VIMS scans along samples from 1 to 35, then increases the line
VIMSexps=[]
for i in range(len(Cass)):
    if Cass[i][5]==Cass[i-1][5] and Cass[i][6]==Cass[i-1][6]+1:
        VIMSexps.append(Cass[i][3]-Cass[i-1][3])
#time to take a pixel in the same line
singlePixExpTime=num.median(VIMSexps)

w=num.where(Cass[:,6]==35)
x=Cass[w][:,3]
#time to swap from one line to the next
lineSwitchTime=num.median((x[1:]-x[:len(x)-1]))-36*singlePixExpTime

#setup the time array. Pick as Tnot the first time in the Cass array
(A,B)=imData.shape
pixelTimes=num.zeros([A,B]).astype(num.float64)

t=0.0
for i in range(A):
    for j in range(B):
        pixelTimes[i,j]=t
        t+=singlePixExpTime
    t+=lineSwitchTime

inot=int(Cass[0,5])
jnot=int(Cass[0,6])
Tnot=Cass[0,3]
pixelTimes+=Tnot-pixelTimes[inot,jnot]


#distance units are in km, resolution units in km/pix
[Xnot,Ynot,Znot]=latLongObj[inot][jnot]['SpacecraftPosition']
distancenot=latLongObj[inot][jnot]['TargetCenterDistance'] 
sampleResolutionnot=latLongObj[inot][jnot]['SampleResolution']/1000.


#velocity polynomial fitting order
polyOrder=2
objx=poly(Cass[:,3],Cass[:,0],renormX=True,order=polyOrder)
objy=poly(Cass[:,3],Cass[:,1],renormX=True,order=polyOrder)
objz=poly(Cass[:,3],Cass[:,2],renormX=True,order=polyOrder)

spaceCraftX=pixelTimes*0.0
spaceCraftY=pixelTimes*0.0
spaceCraftZ=pixelTimes*0.0

for i in range(len(pixelTimes)):
    spaceCraftX[i,:]=objx.eval(pixelTimes[i,:])
    spaceCraftY[i,:]=objy.eval(pixelTimes[i,:])
    spaceCraftZ[i,:]=objz.eval(pixelTimes[i,:])
spaceCraftVectornot=num.array([Xnot,Ynot,Znot])
deltaX=spaceCraftX-spaceCraftX[inot,jnot]
deltaY=spaceCraftY-spaceCraftY[inot,jnot]
deltaZ=spaceCraftZ-spaceCraftZ[inot,jnot]

long_o_not=latLongObj[inot][jnot]['SubSpacecraftLongitude']
lat_o_not=latLongObj[inot][jnot]['SubSpacecraftLatitude']
az_o_not=latLongObj[inot][jnot]['SpacecraftAzimuth']
long_s=latLongObj[inot][jnot]['SubSolarLongitude']
lat_s=latLongObj[inot][jnot]['SubSolarLatitude']
az_s=latLongObj[inot][jnot]['SubSolarAzimuth']

offsets=num.array([130.,140.]) #shift in vertical, shift in horizontal
offsetVel=0.00 #km/s
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

if doFits:
    x,y=0.,0.
    
    nDim=9
    nWalkers=nDim*4
    nBurn=100
    nStep=100

    steps=0
    offWidth=20.
    angWidth=10.
    velWidth=0.1
    distanceWidth=distancenot*0.1
    resWidth=sampleResolutionnot*0.1
    
    r0=[]
    widthArr=num.array([angWidth,angWidth,angWidth,angWidth,angWidth,angWidth,
                        distanceWidth,offWidth,offWidth])
    for i in range(nWalkers):
        entry=num.array([long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,
                         distancenot,offsets[0],offsets[1]])
        entry+=sci.randn(nDim)*widthArr
        #entry[10]=entry[10]%360. #making sure the offset velocity angle is between 0 and 360
        r0.append(entry)
    r0=num.array(r0)
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, callShapeGen, args=[vertices,vertIndices,offsetVel,offsetVelAngle,sampleResolutionnot,inot,jnot,deltaX,deltaY,deltaZ,lons,lats,pixelTimes,imData])

    pos, prob, state = sampler.run_mcmc(r0, nBurn)

    sampler.reset()

    pos, prob, state = sampler.run_mcmc(pos, nStep, rstate0=state)

    samps=sampler.chain
    probs=sampler.lnprobability
    cubeName='junk'
    with open('%s.fit_pickle'%(cubeName),'w+') as outHan:
        pick.dump([samps,probs],outHan)
    
    print sampler.chain
    print
    print sampler.lnprobability

    (Y,X)=probs.shape
    goodSamps=[]
    for i in range(Y):
        for j in range(int(X*0.5),X):
            xx=[]
            for k in range(nDim):
                xx.append(samps[i,j][k])
            xx.append(probs[i,j])
            goodSamps.append(xx[:])
            #goodSamps.append([samps[i,j][0],
            #                  samps[i,j][1],
            #                  samps[i,j][2],
            #                  samps[i,j][3],
            #                  samps[i,j][4],
            #                  samps[i,j][5],
            #                  samps[i,j][6],
            #                  samps[i,j][7],
            #                  samps[i,j][8],
            #                  probs[i,j]])
    goodSamps=num.array(goodSamps)
    (l,ll)=goodSamps.shape
    args=num.argsort(goodSamps[:,ll-1])
    goodSamps=goodSamps[args]

    
    print l,ll
    print 'Best Point: ',goodSamps[l-1]
    print callShapeGen(goodSamps[l-1][:ll-1],vertices,vertIndices,offsetVel,offsetVelAngle,sampleResolutionnot,inot,jnot,deltaX,deltaY,deltaZ,lons,lats,pixelTimes,imData)
    
    sys.exit()
    





###modelling free parameters are:
# space craft distance
# subspacecraft long,lat,az
# subsolar long, lat, az
# y,z image offsets in km
# offsetVel in km/s
# offsetVelAngle in degrees
###Chosen parameters are:
# polyOrder for the velocity modelling
#IR
[long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offX,offY,chi]=[  3.56165447e+02,  -4.15414568e+01,   4.62712643e+01,   7.62419760e+01,
  -1.15523618e+01,   1.70399539e+02,   1.81620185e+04,   1.14675932e+02,
   1.40134996e+02,  -41.]
#Vis
[long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offX,offY,chi]=[  3.50244948e+02,   1.61861891e+00,   4.67986885e+01,   8.00298375e+01,
  -2.14457114e+01,   8.61919423e+01,   1.02140002e+04,   9.76095180e+01,
   2.65240659e+02,  -38]

offsets=num.array([offX,offY])
(image,poly3d,colours,rot_vert,chi)=shapeGen_VIMS(vertices,vertIndices,
                                      long_o_not,lat_o_not,az_o_not,
                                      long_s,lat_s,az_s,
                                      distancenot,
                                      offsets, offsetVel,offsetVelAngle,
                                      sampleResolutionnot,
                                      inot,jnot,
                                      deltaX,deltaY,deltaZ,
                                      lons,lats,
                                      pixelTimes,
                                      imData)
steps=0

callShapeGen((long_o_not,lat_o_not,az_o_not,long_s,lat_s,az_s,distancenot,offsets[0],offsets[1]),vertices,vertIndices,offsetVel,offsetVelAngle,sampleResolutionnot,inot,jnot,deltaX,deltaY,deltaZ,lons,lats,pixelTimes,imData)














 

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


realImPlot=ax3.imshow(256.-imData,aspect='equal',interpolation='nearest')
realImPlot.set_cmap('Greys')

modelImPlot=ax4.imshow(256.-image*256.,aspect='equal',interpolation='nearest')
modelImPlot.set_cmap('Greys')


pyl.show()

sys.exit()
             


#plot!
fig=pyl.figure(1,figsize=(10,10))
fig.subplots_adjust(wspace=0,hspace=0)
ax1=fig.add_subplot(2,2,1,projection='3d')
ax2=fig.add_subplot(2,2,2,projection='3d')
ax3=fig.add_subplot(2,2,3)
ax4=fig.add_subplot(2,2,4)
fig.canvas.mpl_connect('key_press_event',event_handler)


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

#azim,elev=0.0,0.0 gives x+ up, z+ right ,y+ out of screen
#azim,elev=90,0.0 gives x+ up, y+ left, z+ into screen
#azim is a rotation around x.
#elevation is a rotation around z.
#default sbmt orientation is y up, z out of screen, x right
azim=0.0
elev=0.0
#azim=elev=0 and switching x/y and swapping negatives gives correct orientation in sbmt with
#    +x (pylab)=+y sbmt
#    +z (pylab)=+x sbmt
#    +y (pylab)=+z sbmt
ax1.view_init(azim=azim, elev=elev)
ax2.view_init(azim=azim, elev=elev)

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')

key_event=backend_bases.KeyEvent('start',fig.canvas,'d')
event_handler(key_event)
pyl.show()

