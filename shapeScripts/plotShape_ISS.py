#! /usr/bin/env python

from shapeLib import shapeGen_ISS,rot,d2r


import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as pyl
from matplotlib import backend_bases
import numpy as num, scipy as sci
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sys,os
import illuminate as ill
import pickle as pick
import random
from astropy.io import fits as pyf
import time
import emcee
from scipy import signal

def parseCamInfo(fn):
    with open(fn) as cfn:
        data=cfn.readlines()

    out={'SubSpacecraftLongitude': None,
         'SubSpacecraftLatitude': None,
         'SubSpacecraftAzimuth': None,
         'SubSolarLongitude':None,
         'SubSolarLatitude': None,
         'SubSolarAzimuth': None,
         'SampleResolution': None}
    for ii in out:
        for jj in range(len(data)):
            if ii in data[jj]:
                s=data[jj].split()
                out[ii]=float(s[2])
    return out.copy()

def callShapeGen(r,sampleResolution,vertices,vertIndices,imData):
    global steps
    (L_o,l_o,A_o,L_s,l_s,a_s,x,y)=r


    steps+=1
    if abs(x)>60. or abs(y)>60. or x==num.inf or x==-num.inf or x==num.inf or y==-num.inf or L_o>360. or L_o<0. or l_o>90. or l_o<-90. or a_o>360. or a_o<0. or l_s>90. or l_s<-90. or L_s>360. or L_s<0. or a_s>360. or a_s<0.:
        
        return -num.inf

    print steps,L_o,l_o,A_o,L_s,l_s,a_s,x,y,
    
    offsets=num.array([x,y])

    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)=shapeGen_ISS(vertices,vertIndices,
                                                                                     long_o=L_o,lat_o=l_o,az_o=A_o,
                                                                                     long_s=L_s,lat_s=l_s,az_s=a_s,
                                                                                     offsets=offsets,sampleResolution=sampleResolution,
                                                                                     imData=imData)
    print chi
    return 0.5*chi


def callShapeGen_offsetsOnly(r,L_o,l_o,A_o,L_s,l_s,a_s,sampleResolution,vertices,vertIndices,imData):
    global steps
    (x,y)=r


    steps+=1
    if abs(x)>60. or abs(y)>60. or x==num.inf or x==-num.inf or x==num.inf or y==-num.inf or L_o>360. or L_o<0. or l_o>90. or l_o<-90. or a_o>360. or a_o<0. or l_s>90. or l_s<-90. or L_s>360. or L_s<0. or a_s>360. or a_s<0.:
        
        return -num.inf

    print steps,L_o,l_o,A_o,L_s,l_s,a_s,x,y,
    
    offsets=num.array([x,y])

    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)=shapeGen_ISS(vertices,vertIndices,
                                                                                     long_o=L_o,lat_o=l_o,az_o=A_o,
                                                                                     long_s=L_s,lat_s=l_s,az_s=a_s,
                                                                                     offsets=offsets,sampleResolution=sampleResolution,
                                                                                     imData=imData)
    print chi
    return 0.5*chi
 


def event_handler(event):
    global vertics,vertIndices,showSurface,showNormals,showSources
    global long_o,lat_o,a_o,long_s,lat_s,a_s
    global sampleResolution
    global offsets
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

    if event.key=='right':  offsets[0]+=1
    elif event.key=='left': offsets[0]-=1
    elif event.key=='up':   offsets[1]+=1
    elif event.key=='down': offsets[1]-=1
        
    print long_o%360,lat_o,a_o
    print long_s%360.,lat_s,a_s
    print offsets
    
    if event.key=='r':
        long_s+=5.*(random.random()-0.5)
        lat_s+=5.*(random.random()-0.5)
        long_o+=5.*(random.random()-0.5)
        lat_o+=5.*(random.random()-0.5)
        a_o+=5.*(random.random()-0.5)
        #offsets+=(num.array([random.random(),random.random(),random.random()])-0.5)*30.
    
        
    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)=shapeGen_ISS(vertices,vertIndices,
                                                                                     long_o=long_o,lat_o=lat_o,az_o=a_o,
                                                                                     long_s=long_s,lat_s=lat_s,az_s=a_s,
                                                                                     offsets=offsets,sampleResolution=sampleResolution,
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
        vertSpacing=1#0
    else:
        vertSpacing=1
    ax2.scatter(rot_vertices[:,0][::vertSpacing],
                rot_vertices[:,1][::vertSpacing],
                rot_vertices[:,2][::vertSpacing],
                marker='.',s=2,c='k',alpha=0.5)

    if showSources:
        sunSpot=rot_n_sun/num.sum(rot_n_sun**2)**0.5
        ax2.scatter(sunSpot[0]*150.,sunSpot[1]*150.,sunSpot[2]*150.,marker=(5,0,0),c='y',s=70,zorder=0)
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

    ax3.set_xlabel('Simulated Image')
    ax4.set_xlabel('VIMS Image')

    pyl.draw()


#image name and directory it is in
#cubeName='cN1465668217_2'
cubeName='cN1465662798_2'
dir='1465650156_1465674412'

if dir not in os.getcwd():
    os.chdir(dir)
    
#create the fits file
comm='isis2fits from=%s.cub to=junk.fits'%(cubeName)
print comm
print os.system(comm)
print

#load the image data
with pyf.open('/Volumes/data/ISS/coiss_2003/procdata/%s/junk.fits'%(dir)) as han:
    imData=han[0].data

W=num.where(imData<-10000)
imData[W]=num.median(imData)
maxi=num.max(imData)
W=num.where(imData<maxi)
maxi=num.max(imData[W])
imData=num.clip(imData,num.min(imData),maxi)


#(A,B)=imData.shape
#(n,b,p)=pyl.hist(imData.reshape(A*B),bins=500)
#pyl.show()
#sys.exit()
w=num.where(imData>0.001)
imData[w]=256.
W=num.where(imData<256.)
imData[W]=0.0

width=300
imData=imData#[512-width/2:512+width/2,512-width/2:512+width/2]




#########plotting parameters
#Plotting angles and settings
caminfo=parseCamInfo('%s.caminfo'%(cubeName))
long_o=caminfo['SubSpacecraftLongitude']
lat_o=caminfo['SubSpacecraftLatitude']
a_o=caminfo['SubSpacecraftAzimuth']
long_s=caminfo['SubSolarLongitude']
lat_s=caminfo['SubSolarLatitude']
a_s=caminfo['SubSolarAzimuth']
sampleResolution=caminfo['SampleResolution']/1000.
offsets=num.array([0.,0.])
print long_o,lat_o,a_o
print long_s,lat_s,a_s
print



##best point cN1465662789
##chi=-12510.0
long_o=85.00995932
lat_o =-13.70312373
a_o   =253.4702181
long_s=176.44709816
lat_s =-22.3388434
a_s   =81.0656944
offsets=num.array([-10.77297219,    -7.59388273])

##best point cN1465668217_2
##chi=-10709
#long_o=2.76977640e+01
#lat_o =-1.38872320e+01
#a_o   =2.53378246e+02
#long_s=1.24574841e+02
#lat_s =-1.66265355e+01
#a_s   =9.38729824e+01
#offsets=num.array([-1.37270582e+01,  -2.45790072e+00])


##########




showSurface=True
showSources=True
showNormals=True
plotUpdate=False
alpha=0.2
doFits=False
doOffsetFits=False


#load the shape model

###use the following as reference for plotting
#http://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
if not doFits:
    x=1
else:
    x=4
n=(64*x)+1
with open('/Volumes/data/PhoebeShapeModel/CO_SA_ISSNA_5_PHOEBESHAPE_V2_0/data/phoebe_ver%sq.tab'%(n-1)) as han:
    data=han.readlines()

vertices=[]
for i in range(1,6*n**2+1):
    s=data[i].split()
    v=[float(s[1]),float(s[2]),float(s[3])]
    vertices.append(v)

vertices=num.array(vertices)


#swap play
#y=vertices[:,1]*1.
#z=vertices[:,2]*1.
#vertices[:,1]=y*-1.
#vertices[:,2]=z*-1.


#get the vertices indices
vertIndices=[]
for i in range(6*n**2+2,len(data)):
    s=data[i].split()
    vertIndices.append([int(float(s[1]))-1,int(float(s[2]))-1,int(float(s[3]))-1])
vertIndices=num.array(vertIndices)





if doOffsetFits:
    x,y=0.,0.
    
    nDim=2
    nWalkers=20#nDim*8
    nBurn=100
    nStep=500

    steps=0
    offWidth=20.
    angWidth=10.
    r0=[]
    widthArr=num.array([offWidth,offWidth])
    for i in range(nWalkers):
        r0.append(num.array([offsets[0],offsets[1]])+sci.randn(nDim)*widthArr)
    r0=num.array(r0)
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, callShapeGen_offsetsOnly, args=[long_o,lat_o,a_o,long_s,lat_s,a_s,sampleResolution,vertices,vertIndices,imData])

    pos, prob, state = sampler.run_mcmc(r0, nBurn)

    sampler.reset()

    pos, prob, state = sampler.run_mcmc(pos, nStep, rstate0=state)

    samps=sampler.chain
    probs=sampler.lnprobability

    with open('%s.offsetfit_pickle'%(cubeName),'w+') as outHan:
        pick.dump([samps,probs],outHan)
    
    print sampler.chain
    print
    print sampler.lnprobability

    (Y,X)=probs.shape
    goodSamps=[]
    for i in range(Y):
        for j in range(int(X*0.5),X):
            goodSamps.append([samps[i,j][0],
                              samps[i,j][1],
                              samps[i,j][2],
                              samps[i,j][3],
                              samps[i,j][4],
                              samps[i,j][5],
                              samps[i,j][6],
                              samps[i,j][7],
                              probs[i,j]])
    goodSamps=num.array(goodSamps)
    args=num.argsort(goodSamps[:,8])
    goodSamps=goodSamps[args]

    
    (l,ll)=goodSamps.shape
    print 'Best Point: ',goodSamps[l-1]
    print callShapeGen(goodSamps[l-1][:8],sampleResolution,vertices,vertIndices,imData)
    
    sys.exit()
    


if doFits:
    x,y=0.,0.
    
    nDim=8
    nWalkers=100#nDim*8
    nBurn=200
    nStep=200

    steps=0
    offWidth=20.
    angWidth=10.
    r0=[]
    widthArr=num.array([angWidth,angWidth,angWidth,angWidth,angWidth,angWidth,offWidth,offWidth])
    for i in range(nWalkers):
        r0.append(num.array([long_o,lat_o,a_o,long_s,lat_s,a_s,offsets[0],offsets[1]])+sci.randn(nDim)*widthArr)
    r0=num.array(r0)
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, callShapeGen, args=[sampleResolution,vertices,vertIndices,imData])

    pos, prob, state = sampler.run_mcmc(r0, nBurn)

    sampler.reset()

    pos, prob, state = sampler.run_mcmc(pos, nStep, rstate0=state)

    samps=sampler.chain
    probs=sampler.lnprobability

    with open('%s.fit_pickle'%(cubeName),'w+') as outHan:
        pick.dump([samps,probs],outHan)
    
    print sampler.chain
    print
    print sampler.lnprobability

    (Y,X)=probs.shape
    goodSamps=[]
    for i in range(Y):
        for j in range(int(X*0.5),X):
            goodSamps.append([samps[i,j][0],
                              samps[i,j][1],
                              samps[i,j][2],
                              samps[i,j][3],
                              samps[i,j][4],
                              samps[i,j][5],
                              samps[i,j][6],
                              samps[i,j][7],
                              probs[i,j]])
    goodSamps=num.array(goodSamps)
    args=num.argsort(goodSamps[:,8])
    goodSamps=goodSamps[args]

    
    (l,ll)=goodSamps.shape
    print 'Best Point: ',goodSamps[l-1]
    print callShapeGen(goodSamps[l-1][:8],sampleResolution,vertices,vertIndices,imData)
    
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

ax1.view_init(azim=azim, elev=elev)
ax2.view_init(azim=azim, elev=elev)

ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.set_zlabel('z')

key_event=backend_bases.KeyEvent('start',fig.canvas,'d')
event_handler(key_event)
pyl.show()

