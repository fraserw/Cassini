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
import pickle as pick
import random
from astropy.io import fits as pyf
import time
import emcee
from scipy import signal


d2r=num.pi/180.
r2d=180./num.pi
twopi=2*num.pi



def rot(ang,axis):
    if axis=='x':
        w=num.array([[1.0,0.0,0.0],
                     [0.0,num.cos(ang),-num.sin(ang)],
                     [0.0,num.sin(ang),num.cos(ang)]])
    elif axis=='y':
        w=num.array([[num.cos(ang),0.0,-num.sin(ang)],
                     [0.0,1.0,0.0],
                     [num.sin(ang),0.0,num.cos(ang)]])
    elif axis=='z':
        w=num.array([[num.cos(ang),-num.sin(ang),0.0],
                     [num.sin(ang),num.cos(ang),0.0],
                     [0.0,0.0,1.0]])
    return w

def arrayMagRep(inVec):
    xxx=inVec.shape
    (AA,BB,CC)=inVec.shape
    x=inVec**2
    mags=num.sum(x,axis=2)**0.5
    repMags=num.repeat(mags,CC,axis=1).reshape(AA,BB,CC)
    return repMags
def arrayMag2(inVec):
    (AA,BB)=inVec.shape
    x=inVec**2
    mags=num.sum(x,axis=1)**0.5
    return mags


def reshape(vertices,vertIndices,
            phi_o=200.35216011092,theta_o=-21.76086293882,az_o=264.28748413648,
            phi_s=291.5572425051,theta_s=-12.85085685023):
    

    #Pass the actual camInfo log values (in degrees and with the sign in the log file!)
    #
    #        
    #defaults are from ISS image cN1465652275_1.caminfo
    ####it seems we need negatives infront of all the angles!
    #longitude phi_o=d2r*-(200.35216011092)
    #latitude theta_o=d2r*-(-21.76086293882)
    #azimuth az_o=d2r*-( 264.28748413648)



    
    Aphi_s=phi_s*-d2r
    Aphi_o=phi_o*-d2r
    Atheta_s=theta_s*-d2r
    Atheta_o=theta_o*-d2r
    Aaz_o=az_o*-d2r

    n_obs_r=num.array([num.cos(Aphi_o)*num.cos(Atheta_o),
                      num.sin(Aphi_o)*num.cos(Atheta_o),
                      num.sin(Atheta_o)])
    n_sun=num.array([num.cos(Aphi_s)*num.cos(Atheta_s),
                    num.sin(Aphi_s)*num.cos(Atheta_s),
                    num.sin(Atheta_s)])*1.e15

    n_obs=num.array([1.,0.,0.])*1.e15


    #apply a triplet of rotations to the body to get the geometry correct
    #longitude pivots around the z axis and has zero at +i (xhat)
    #latitude pivots around the y axis and has zero in the xy plane
    #azimuth is about the z axis again

    w1=rot(Aphi_o,'z')
    w2=rot(Atheta_o,'y')
    w3=rot(Aaz_o,'x')
    W=num.dot(w3,num.dot(w2,w1))

    w1=rot(-Aphi_o,'z')
    w2=rot(-Atheta_o,'y')
    nW=num.dot(w2,w1)

    
    rot_n_sun=num.dot(nW,n_sun)

     
    #rot_vertices=vertices*0.0
    #for i in range(len(vertices)):
    #    rot_vertices[i]=num.dot(W,vertices[i])
    rot_vertices=[]
    for i in range(len(vertices)):
        rot_vertices.append(num.dot(W,vertices[i]))
    rot_vertices=num.array(rot_vertices)

        
    #convert the shape facets to the polygon needed by matplotlib
    #vertIndices=[]
    #for i in range(6*65**2+2,len(data)):
    #    s=data[i].split()
    #    vertIndices.append([int(float(s[1]))-1,int(float(s[2]))-1,int(float(s[3]))-1])
    #vertIndices=num.array(vertIndices)
    poly3d=num.copy(rot_vertices[vertIndices])



    ####making getNormalMidpoint faster!
    ####to the next #### replaces illuminate.getNormalMidpoint
    rotMids=num.mean(poly3d,axis=1) #midpoints

    (A,B)=rotMids.shape
    repMids=num.repeat(rotMids,3,axis=0).reshape(A,B,3) 
    abc=poly3d-repMids  #non-normalized abc vectors
    abc/=arrayMagRep(abc)  #abc now normzalized

    a=abc[:,0]
    b=abc[:,1]
    c=abc[:,2]

    axb=num.cross(a,b)
    bxc=num.cross(b,c)
    cxa=num.cross(c,a)

    cross=num.concatenate([num.array([axb]).transpose(),num.array([bxc]).transpose(),num.array([cxa]).transpose()]).transpose().reshape(A,3,3) #unnormalized cross product vectors
    cross/=arrayMagRep(cross) #normalized cross products
    
    normals=num.mean(cross,axis=1) #normals

    #print rotMids,normals.shape
    w=num.where(arrayMag2(rotMids-normals)>arrayMag2(rotMids+normals))
    normals[w]*=-1

    n=normals
    m=rotMids
    #### end replacement
    
    #get the observability
    #illuminated, visible, and both facets
    mag_n_sun=num.sum(rot_n_sun**2)**0.5
    mag_n_obs=num.sum(n_obs**2)**0.5
    
    angs_s=num.arccos(num.dot(n,rot_n_sun)/mag_n_sun)
    illuminated=num.where((angs_s<=num.pi/2)&(angs_s>=-num.pi/2))[0]

    angs_o=num.arccos(num.dot(n,n_obs)/mag_n_obs)
    visible=num.where((angs_o>=-num.pi/2)&(angs_o<=num.pi/2))[0]

    ill_and_obs=num.where( ((angs_s>=-num.pi/2)&(angs_s<=num.pi/2))
                        & ((angs_o>=-num.pi/2)&(angs_o<=num.pi/2)))[0]

    not_ill_and_obs=num.where(((angs_s>=num.pi/2)|(angs_s<=-num.pi/2))
                            &((angs_o>=-num.pi/2)&(angs_o<=num.pi/2)))[0]



    #setup the colours array
    #0.0 is black,1.0 is white
    ###since the fitting code doesn't want colour but just illuminated or not, we could just do colours=1 where ill_and_obs==True, may be a little faster
    colours=num.repeat(num.cos(angs_o)*num.cos(angs_s),3).reshape(len(m),3)
    colours[num.where(colours<0)]=0.0
    colours[not_ill_and_obs]=0.0
    cmax=num.max(colours)
    colours/=cmax



    #collection=Poly3DCollection(poly3d[ill_and_obs],linewidths=0.0,facecolors=colours[ill_and_obs])
    #collection.set_alpha(alpha)
    poly2d=num.copy(poly3d[ill_and_obs])
    colours2d=num.copy(colours[ill_and_obs])

    l1=((poly2d[:,0,0]-poly2d[:,1,0])**2+(poly2d[:,0,1]-poly2d[:,1,1])**2)**0.5
    l2=((poly2d[:,0,0]-poly2d[:,2,0])**2+(poly2d[:,0,1]-poly2d[:,2,1])**2)**0.5
    l3=((poly2d[:,1,0]-poly2d[:,2,0])**2+(poly2d[:,1,1]-poly2d[:,2,1])**2)**0.5
    s=(l1+l2+l3)/2.
    area=(s*(s-l1)*(s-l2)*(s-l3))**0.5
    area[num.where(num.isnan(area))]=0.


    return (rot_vertices,num.copy(poly3d[ill_and_obs]),num.copy(colours[ill_and_obs]),num.copy(m[ill_and_obs]),colours2d,area,n,m,n_obs,rot_n_sun,angs_s,angs_o)


def callShapeGen(r,vertices,vertIndices,imData):
    global steps
    (L_o,l_o,A_o,L_s,l_s,x,y)=r

    A_s=0.0
    steps+=1
    if abs(x)>60. or abs(y)>60. or x==num.inf or x==-num.inf or x==num.inf or y==-num.inf or L_o>360. or L_o<0. or l_o>90. or l_o<-90. or a_o>360. or a_o<0. or l_s>90. or l_s<-90. or L_s>360. or L_s<0.:
        
        return -num.inf

    print steps,L_o,l_o,A_o,L_s,l_s,x,y,
    

    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)=shapeGen(vertices,vertIndices,
                                                                                     long_o=L_o,lat_o=l_o,az_o=A_o,
                                                                                     long_s=L_s,lat_s=l_s,
                                                                                     imData=imData)
    print chi
    return 0.5*chi
 

def shapeGen(vertices,vertIndices,
             long_o,lat_o,az_o,
             long_s,lat_s,
			 planetVertex,
             imData):

    x=reshape(vertices,vertIndices,
                phi_o=long_o,theta_o=lat_o,az_o=az_o,
                phi_s=long_s,theta_s=lat_s)

    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o)=x

    rotCentreSpot=m2d[num.where(m2d==m[planetVertex])]

    print
    print m,rotCentreSpot,m[planetVertex]
    print
    w=num.where((num.abs(m2d[:,1]-rotCentreSpot[1])<(sampleResolution/2.)) & (num.abs(m2d[:,2]-rotCentreSpot[2])<(sampleResolution/2.)))

    print num.sum(c2d[w])
    return
    sys.exit()
    #print m2d[w]
    #print rotCentreSpot
    #print sampleResolution


    

    #generate the model image
    c2d=c2d[:,0]
    c2d*=a2d
    c2d/=num.max(c2d)
    print c2d
    sys.exit()

    (A,B)=imData.shape
    y=num.arange(-A*sampleResolution/2,A*sampleResolution/2,sampleResolution)
    z=num.arange(-B*sampleResolution/2,B*sampleResolution/2,sampleResolution)
    image=num.zeros([len(y),len(z)]).astype('float')
    
    K=((m2d[:,1]-y[0])/(y[1]-y[0])).astype('int')
    L=((m2d[:,2]-z[0])/(z[1]-z[0])).astype('int')

    
    w=num.where((K>=0)&(K<A)&(L>=0)&(L<B))

    Kw=K[w]
    Lw=L[w]
    c2dw=c2d[w]


    for ii in range(len(Kw)):
        image[Kw[ii],Lw[ii]]+=c2dw[ii]

    imMax=num.max(image)
    #debugging if len(w[0])>0: print '\n',imMax,num.max(c2dw),'***'
    image/=imMax
    image*=256.
    image=image.transpose()[::-1,:]

    
    w=num.where(image>0)
    image[w]=256.
    #done image generation
    

    #256 is illuminated, 0 is dark
    chi=-num.sum(((imData/256.)-(image/256.))**2)#/(A*B-7.)
    #print chi


    return (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)

            
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
    
        
    (rot_vertices,p3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o,chi,image)=shapeGen(vertices,vertIndices,
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
    
image='2004163T121836_2004163T192848/cv1465671822_1_ir'

with pyf.open('/Users/fraserw/data/VIMS/covims_0004/procdata/%s_mean.fits'%(image)) as han:
    imData=han[0].data

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

width=300
imData=imData

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

mids=num.mean(vertices[vertIndices],axis=1) #midpoints
lons=num.arctan2(mids[:,1],mids[:,0])*r2d
lats=num.arcsin(mids[:,2]/(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5)*r2d


#line is i, sample is j
for i in range(len(latLongObj)): 
    for j in range(len(latLongObj[i])):
        if latLongObj[i][j]['SubSolarLongitude']<>None:

            #########plotting parameters
            #Plotting angles and settings
            print i+1,j+1,
            long_o=latLongObj[i][j]['SubSpacecraftLongitude']
            lat_o=latLongObj[i][j]['SubSpacecraftLatitude']
            az_o=latLongObj[i][j]['SpacecraftAzimuth']-180.
            long_s=latLongObj[i][j]['SubSolarLongitude']
            lat_s=latLongObj[i][j]['SubSolarLatitude']
    
            sampleResolution=latLongObj[i][j]['SampleResolution']/1000.

            #print long_o,lat_o,az_o
            #print long_s,lat_s
            #print sampleResolution
            #sys.exit()    
            planetLat=latLongObj[i][j]['PlanetocentricLatitude']
            eastLong=latLongObj[i][j]['PositiveEast360Longitude']
            ##########
    
            planetVertex=num.argmin(((planetLat-lats)**2+(eastLong-lons)**2)**0.5)

            print planetLat,eastLong
            shapeGen(vertices,vertIndices,
                    long_o,lat_o,az_o,
                    long_s,lat_s,
                    planetVertex,
                    imData)


if doFits:
    x,y=0.,0.
    
    nDim=7
    nWalkers=48
    nBurn=200
    nStep=100

    steps=0
    offWidth=15.
    angWidth=20.
    r0=[]
    widthArr=num.array([angWidth,angWidth,angWidth,angWidth,angWidth,offWidth,offWidth])
    for i in range(nWalkers):
        r0.append(num.array([long_o,lat_o,a_o,long_s,lat_s,offsets[0],offsets[1]])+sci.randn(nDim)*widthArr)
    r0=num.array(r0)
    
    sampler = emcee.EnsembleSampler(nWalkers, nDim, callShapeGen, args=[vertices,vertIndices,imData])

    pos, prob, state = sampler.run_mcmc(r0, nBurn)

    sampler.reset()

    pos, prob, state = sampler.run_mcmc(pos, nStep, rstate0=state)

    samps=sampler.chain
    probs=sampler.lnprobability
    
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
                              probs[i,j]])
    goodSamps=num.array(goodSamps)
    args=num.argsort(goodSamps[:,7])
    goodSamps=goodSamps[args]

    
    (l,ll)=goodSamps.shape
    print 'Best Point: ',goodSamps[l-1]
    print callShapeGen(goodSamps[l-1][:7],vertices,vertIndices,imData)
    print 'Best Point: ',goodSamps[0]
    print callShapeGen(goodSamps[0][:7],vertices,vertIndices,imData)
    
    print num.max(probs[:,int(X*0.5):]),
    print num.min(probs[:,int(X*0.5):]),'**'

    sys.exit()
    
#(long_o,lat_o,az_o,lat_s,az_s,x,y,z)=r



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

