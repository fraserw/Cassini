#! /usr/bin/env python

from shapeLib import shapeGen_ISS,rot,d2r,midsNormals,arrayMag2,r2d
import sys

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

def closeWithg(event):
    if event.key=='g': pyl.close()
        
def goodCrater(event):
    global good
    if event.key in ['b','B']:
        good=False
        print 'Not a crater,\n'
    pyl.close()
    
def getRadius(event):
    global selRad
    if event.button==1:
        print "Selected radius %s."%(event.xdata)
        selRad=event.xdata
    elif event.button==3:
        print 'Marked as not a crater!'
    pyl.close()

def event_handler(event):
    global vertics,vertIndices,showSurface,showNormals,showSources
    global long_o,lat_o,a_o,long_s,lat_s,a_s
    global sampleResolution
    global offsets
    global plotUpdate
    global alpha
    global azim,elev
    global www,Parg

    print event.key

    imData=num.ones([10,10])
    
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
        #determine under point
        w=num.where(m[:,0]>0)[0]
        arg=num.argmin(m[w][:,1]**2+m[w][:,2]**2)
        w=num.where((m[w[arg]][0]==m[:,0])&(m[w[arg]][1]==m[:,1])&(m[w[arg]][2]==m[:,2]))[0]

        under_arg=w[0]

        #print m[under_arg]
        distFromUnder=((m[:,0]-m[under_arg][0])**2 + (m[:,1]-m[under_arg][1])**2 +(m[:,2]-m[under_arg][2])**2)**0.5
        w=num.where(distFromUnder<8.)
        #print m[0],m[1]
        #print w
        print www
        
        
        for i in www:#0,len(n),45):
            #if angs_s[i]<=num.pi/2.:
            Y=num.concatenate([m[i],m[i]+n[i]*50.]).reshape(2,3)
            ax2.plot(Y[:,0],Y[:,1],Y[:,2],c='r',lw=3)
            if i==Parg:            ax2.plot(Y[:,0],Y[:,1],Y[:,2],c='b',lw=3)

    ax1.view_init(azim=azim, elev=elev)
    ax2.view_init(azim=azim, elev=elev)
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')

    ax3.set_xlabel('Simulated Image')
    ax4.set_xlabel('VIMS Image')

    pyl.draw()



#########plotting parameters
#Plotting angles and settings
long_o=112.8
lat_o=14.6
a_o=0.0
long_s=190.0
lat_s=14.6
a_s=0.0
sampleResolution=1.
offsets=num.array([0.,0.])



showSurface=True
showSources=True
showNormals=True
plotUpdate=False
alpha=0.2


#load the shape model

###use the following as reference for plotting
#http://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
#shape model resolution
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


#get the vertices indices
vertIndices=[]
for i in range(6*n**2+2,len(data)):
    s=data[i].split()
    vertIndices.append([int(float(s[1]))-1,int(float(s[2]))-1,int(float(s[3]))-1])
vertIndices=num.array(vertIndices)

(mids,normals)=midsNormals(vertices[vertIndices])
lons=(num.arctan2(mids[:,1],mids[:,0])*r2d)%360
lats=num.arcsin(mids[:,2]/(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5)*r2d

Parg=num.argmin(((lons-190.)**2+(lats-8.)**2)**0.5)


searchRadius=2.
colours=['y','b','g','m','r']
pL=[[],[],[],[],[]]
pl=[[],[],[],[],[]]
taken=[]
craters=[]
fracCurves=[]
cratMids=[]

job = int(float(sys.argv[1]))
for i in range(job*100000,min((job+1)*100000,len(mids))):
    if i in taken: continue
    #print i,mids[i],lons[i],lats[i]
    d=((mids[:,0]-mids[i][0])**2 + (mids[:,1]-mids[i][1])**2 + (mids[:,2]-mids[i,2])**2)**0.5
    www=num.where(d<searchRadius)[0]
    
    rs=mids[www]-mids[i]
    mag_rs=num.repeat(arrayMag2(rs),3).reshape(len(rs),3)
    rs/=mag_rs

    normAngles=[]
    for k in range(len(rs)):
        a=num.arccos(num.dot(normals[www][k],rs[k]))
        if not num.isnan(a):
            normAngles.append(a)
    normAngles=num.array(normAngles)

    w=num.where(normAngles<=num.pi/2.)
    frac=float(len(w[0]))/len(normAngles)
    measureCrater=False
    if frac<0.25:
        measureCrater=True
        #print i,len(mids),mids[i],lons[i],lats[i],frac
        #for k in range(len(www)):
        #    pL[0].append(lons[www[k]])
        #    pl[0].append(lats[www[k]])

    if not measureCrater: continue
    
    inCratFracs=[]
    medPoints=[]
    medarg=i
    cratRadss=num.concatenate([num.arange(2.,4.,0.2),num.arange(4.,8.,0.4),num.arange(8.,15.,0.7),num.arange(15.,25.,2.5)])
    for j in range(len(cratRadss)):
        www=num.where(d<cratRadss[j])[0]
        rs=mids[www]-mids[medarg]
        mag_rs=num.repeat(arrayMag2(rs),3).reshape(len(rs),3)
        rs/=mag_rs
        
        normAngles=[]
        for k in range(len(rs)):
            a=num.arccos(num.dot(normals[www][k],rs[k]))
            if not num.isnan(a):
                normAngles.append(a)
        normAngles=num.array(normAngles)
        w=num.where(normAngles>num.pi/2.)
        medLong=num.median(lons[www[w]])
        medLat=num.median(lats[www[w]])
        medarg=num.argmin(((lons-medLong)**2+(lats-medLat)**2)**0.5)
        medPoints.append([medarg,medLong,medLat])
        d=((mids[:,0]-mids[medarg][0])**2 + (mids[:,1]-mids[medarg][1])**2 + (mids[:,2]-mids[medarg,2])**2)**0.5
        #print cratRadss[j],medarg,medLong,medLat,len(w[0]),len(normAngles),float(len(w[0]))/len(normAngles)
        inCratFracs.append(float(len(w[0]))/len(normAngles))
    medPoints=num.array(medPoints)
    inCratFracs=num.array(inCratFracs)



    arg=num.argmin(inCratFracs)
    argmax=num.argmax(inCratFracs)
    for j in range(arg,0,-1):
        if inCratFracs[j-1]<inCratFracs[j]:
            #print cratRadss[j+1]
            if inCratFracs[j]>=inCratFracs[0]:
                selRad=cratRadss[j+1]
            else:
                selRad=-1 #reject the crater
            break
    if argmax>arg and argmax<>len(cratRadss)-1:
        selRad=cratRadss[argmax+1]
        
    if selRad==-1: continue

    """
    fig=pyl.figure('Crater')
    fig.subplots_adjust(hspace=0)
    sp1=fig.add_subplot(311,xticklabels='')
    pyl.plot([selRad,selRad],[num.min(inCratFracs),num.max(inCratFracs)],'k-',lw=2)
    pyl.plot(cratRadss,inCratFracs)
    sp2=fig.add_subplot(312,xticklabels='')
    pyl.plot(cratRadss,medPoints[:,1])
    sp3=fig.add_subplot(313)
    pyl.plot(cratRadss,medPoints[:,2])
    pyl.connect('key_press_event',closeWithg)
    #pyl.connect('button_press_event',getRadius)
    pyl.show()
    pyl.close()
    """
        
    
    K=num.sum(num.less_equal(cratRadss,selRad))-1
    if K==len(cratRadss): K-=1
    bestRad=cratRadss[K]
    medarg=medPoints[K][0]
    print i,len(mids),medPoints[K]
    if medarg in taken:
        print 'Point already taken'
        continue
    d=((mids[:,0]-mids[medarg][0])**2 + (mids[:,1]-mids[medarg][1])**2 + (mids[:,2]-mids[medarg,2])**2)**0.5
    
    #now identify all the points with the right angle to the mid normal, and trace it
    www=num.where(d<bestRad)[0] 
    rs=mids[www]-mids[medarg]
    mag_rs=num.repeat(arrayMag2(rs),3).reshape(len(rs),3)
    rs/=mag_rs

    normAngles=[]
    for k in range(len(rs)):
        a=num.arccos(num.dot(normals[www][k],rs[k]))
        if not num.isnan(a):
            normAngles.append(a)
    normAngles=num.array(normAngles)
    w=num.where(normAngles>num.pi/2.) #points inside crater
    
    cratMids.append(mids[www[w]])
    #pyl.plot(midsInCrat[:,0],midsInCrat[:,1])
    #pyl.show()
    #pyl.plot(midsInCrat[:,0],midsInCrat[:,2])
    #pyl.show()
    #pyl.plot(midsInCrat[:,1],midsInCrat[:,2])
    #pyl.show()

    good=True
    #pyl.title(medPoints[K])
    #pyl.scatter(lons[www[w]],lats[www[w]],c='r',zorder=10)
    #w=num.where(normAngles<=num.pi/2.)
    #pyl.scatter(lons[www[w]],lats[www[w]],zorder=10)
    #pyl.connect('key_press_event',goodCrater)
    #pyl.show()
    #pyl.close()
    if not good: continue
    
    for j in range(len(www)):
        taken.append(www[j])
   
    craters.append([medPoints[K][0],medPoints[K][1],medPoints[K][2],bestRad])
    fracCurves.append(inCratFracs[:])
    print
    print i,craters
    print

    fn='craters/craters_x%s_%s.pickle'%(str(x),job)
    with open(fn,'w+') as outhan:
        pick.dump([craters,fracCurves,cratMids,taken],outhan)
        
for j in range(len(pL)):
    pyl.scatter(pL[j],pl[j],c=colours[j],zorder=10-j)
pyl.show()
sys.exit()


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

