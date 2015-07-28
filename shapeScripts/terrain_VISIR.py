#! /usr/bin/env pythonw

import numpy as num, pylab as pyl, pickle as pick,os
from astropy.io import fits as pyf
from pixTimes import getPixTimes,getVels
from shapeLib import *
import pixelLatLong as pll
from MChandler import getFit
import specAnalysis
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib.collections import PolyCollection
from matplotlib import colors
import pickle as pick
from pds.imageextractor import ImageExtractor
from pds.core.common import open_pds
from matplotlib import image as mpimg
#from scipy.interpolate import griddata
from scipy import interpolate as interp
import notoriousTEX
from scipy import stats




viswave,irwave=specAnalysis.viswave*1.,specAnalysis.irwave*1.
notoriousTEX.plotSettings(8,8)

class line:
    def __init__(x0,y0,x1,y1):
        self.x0=x0
        self.y0=y0
        self.x1=x1
        self.y1=y1
        self.m=(y1-y0)/(x1-x0)
        self.b=y1-m*x1
    def __call__(self,x):
        return x*self.m+self.b


            
class profileHandler:

    def __init__(self,figure,sp,lons,lats,elevs,mids,water,mixColours,width=1,tol=10):
        """Width is the distance in degrees away from the line before a vertice is considered for plotting."""
        
        self.tol=tol
        self.figure=figure
        self.main_figure_title=self.figure.canvas.get_window_title()
        self.sp=sp
        self.width=width

        self.lons=lons*1.
        self.lats=lats*1.
        self.elev=D*1.
        self.mids=mids
        
        self.water=water*1.
        self.mixColours=mixColours
        
        self.figure.canvas.mpl_connect('pick_event',self.pick)
        self.figure.canvas.mpl_connect('button_release_event',self.release)
        self.spots=num.array([[-17.,20.],[26.,25.]])

        (dx,dy)=self.offset()
         
        self.picked=None

        self.profileFigure=pyl.figure('Profile Figure')
        self.profileAxis=self.profileFigure.add_subplot(111)

        (inside,lengthAlongCut)=self.getInside()
        self._drawProfile(inside,lengthAlongCut)

        
        pyl.figure(self.main_figure_title)
        pyl.draw()

    def offset(self):
        angle=num.arctan2( (self.spots[1][1]-self.spots[0][1]),(self.spots[1][0]-self.spots[0][0]))
        dy=-num.cos(angle)*self.width/2.
        dx=num.sin(angle)*self.width/2.

        try:
            self.lineArtist0.remove()
            self.lineArtist1.remove()
            self.lineArtist2.remove()
        except: pass
            
        self.lineArtist0,=pyl.plot(self.spots[:,0],self.spots[:,1],'k-o',lw=0,picker=self.tol,zorder=10)
        self.lineArtist1,=pyl.plot(self.spots[:,0]-dx,self.spots[:,1]-dy,'k-',lw=2,picker=self.tol,zorder=10)
        self.lineArtist2,=pyl.plot(self.spots[:,0]+dx,self.spots[:,1]+dy,'k-',lw=2,picker=self.tol,zorder=10)
        return (dx,dy)

    def getInside(self):
        A=self.spots[1]-self.spots[0]
        dists=((self.spots[0][0]-lons)**2+(self.spots[0][1]-lats)**2)**0.5
        closestMid=num.argmin(dists)

        lengthAlongCut=num.zeros(len(self.lons))
        inside=num.zeros(len(self.lons)).astype('int')
        for ii in range(len(self.lons)):
            (l,bool)=self._getInside(A,self.lons[ii],self.lats[ii])
            inside[ii]=bool
            lengthAlongCut[ii]=l#*num.cos(self.lats[ii])

            """
            #uses the haversine formula
            dL=self.lons[ii]-self.lons[closestMid]
            dl=self.lats[ii]-self.lats[closestMid]
            ah=num.sin(dl*d2r/2.)**2+num.cos(self.lats[ii]*d2r)*num.cos(self.lats[closestMid])*num.sin(dL*d2r/2.)**2
            angle=2*num.arctan2(ah**0.5,(1.-ah)**0.5)
            arc=106.6*angle
            if bool:print dL,dl,ah,angle,arc
            lengthAlongCut[ii]=arc
            """
            
            """
            if bool:
                a2=num.sum((mids[closestMid])**2)
                b2=num.sum((mids[ii])**2)
                l2=num.sum((mids[closestMid]-mids[ii])**2)
                if lons[ii]/lons[closestMid]<0:
                    angle=num.arccos((l2-a2-b2)/(-2*a2**0.5*b2**0.5))
                else:
                    angle=num.arccos((l2-a2-b2)/(-2*a2**0.5*b2**0.5))
                arc=angle*106.6
                #print a2,b2,l2,angle
                lengthAlongCut[ii]=arc
            """
        w=num.where(inside)[0]
        return (w,lengthAlongCut)
            
    def _getInside(self,A,test_long=1.,test_lat=21.):
        #A=self.spots[1]-self.spots[0]
        B=num.array([test_long,test_lat])-self.spots[0] #the coordinate inside the line?
        magA=(A[0]**2+A[1]**2)**0.5
        magB=(B[0]**2+B[1]**2)**0.5

        a=A*num.dot(A,B)/(magA*magA)
        b=B-a

        maga=(a[0]**2+a[1]**2)**0.5
        magb=(b[0]**2+b[1]**2)**0.5
        
        if magb<self.width/2. and a[0]>=0 and a[0]<=A[0]: return (maga,True)
        else: return (maga,False)


    def _drawProfile(self,inside,lengthAlongCut):
        
        pyl.figure('Profile Figure')
        try:
            self.elev_scat.remove()
        except:
            pass
        
        www=num.where(self.water[inside][:,3]>0)
        self.elev_scat=pyl.scatter(lengthAlongCut[inside][www],self.elev[inside][www],s=60,c=self.water[inside][www],alpha=1.)
        self.profileAxis.set_xlim([num.min(lengthAlongCut[inside]),num.max(lengthAlongCut[inside])])
        self.profileAxis.set_ylim([num.min(self.elev[inside]),num.max(self.elev[inside])])
        pyl.draw()

        
    def pick(self,event):
        (x,y)=(event.mouseevent.xdata,event.mouseevent.ydata)
        d=((x-self.spots[:,0])**2+(self.spots[:,1]-y)**2)**0.5
        arg=num.argmin(d)
        if d[arg]<self.tol:
            self.picked=arg
            self.xy=[x,y]

    def release(self,event):
        x,y=event.xdata,event.ydata
        if self.picked<>None:
            self.spots[self.picked][0]=x
            self.spots[self.picked][1]=y
            self.picked=None
            (self.dx,self.dy)=self.offset()

            (inside,lengthAlongCut)=self.getInside()
            self._drawProfile(inside,lengthAlongCut)

            pyl.figure(self.main_figure_title)
            pyl.draw()


def isValidSpecValues(r):
    if (r[3]<0.0 or r[3]>0.5) or (r[4]<0.55 or r[4]>0.85) or r[3]==None or r[4]==None:
        return False
    else: return True





###for removing trailing zeros on contours
class nf(float):
     def __repr__(self):
         str = '%.1f' % (self.__float__(),)
         if str[-1]=='0':
             return '%.0f' % self.__float__()
         else:
             return '%.1f' % self.__float__()

def normVisSpec(spec,w0=70,w1=86):
    #global viswave,irwave
    #w=num.where((viswave>irwave[0]) & (viswave<0.98))
    #print w
    return num.median(spec[w0:w1+1])
    return spec/num.median(spec[w0:w1+1])
def normIRSpec(spec,w0=0,w1=7):
    #global viswave,irwave
    #w=num.where(irwave<viswave[86])
    #print w
    #sys.exit()
    return num.median(spec[w0:w1+1])
    return spec/num.median(spec[w0:w1+1])

def nanmean(x):
    o=[]
    if x.shape[1]>3:
        for ii in range(len(x[0])):
            w=num.where(num.isnan(x[:,ii])==False)
            o.append(num.mean(x[:,ii][w]))
    else:
        for ii in range(len(x[0])):
            w=num.where(num.isnan(x[:,ii])==False)
            o.append(num.mean(x[:,ii][w]))
    return num.array(o)
            
######
#test comparisons between the PDS image and the coordinate system generated by my lons/lats scripts below shows that they agree.
#used junk.py to generate a crater map, compared that with PhoebeFull.tiff. Good example, in my map, a crater at -116 to -103, which
#corresponds to ~118-103 left of 0, and 244 to 257 in the sbmt display
#all three coordinate systems agree!*******
######
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


def zeroBin(image,binWidth=25):
    (A,B)=image.shape
    mini=0.02#num.min(image)
    outImage=num.zeros([A/binWidth,B/binWidth]).astype('float64')
    for ii in range(0,A-binWidth,binWidth):
        for jj in range(0,B-binWidth,binWidth):
            subsec=image[ii:ii+binWidth,jj:jj+binWidth]
            w=num.where(subsec>mini)
            outImage[ii/binWidth,jj/binWidth]=num.mean(subsec[w])
    return outImage

def getAlbedo(lng,lat,binIm):
    (height,width)=binIm.shape
    x=(  (lng+180)*width/360.).astype('int')
    y=(((lat+90)*height/180.)).astype('int')

    return binIm[y,x]


pyl.rcParams['contour.negative_linestyle'] = 'dashed'


"""
#for the PDS maps
ie = ImageExtractor()
(img, labels) = ie.extract(open_pds('/data/PhoebePDSMaps/SP_1M_0_0_MERCATOR.IMG'))
img.save('Mercator.png')
sys.exit()
"""
newImage=mpimg.imread('/data/PhoebePDSMaps/PhoebeFull.png')
newImage=num.sum(newImage[:,:,0:3],axis=2)/3.
binnedNewImage=zeroBin(newImage)

with pyf.open('/data/PhoebeAlbedoMap/Sim99mod2.fits') as ddd:
    albedoImage=ddd[0].data*1.
w=num.where(albedoImage<115)
W=num.where(albedoImage>=115)
albedoImage-=115.
albedoImage*=0.07/(221.-115.)
albedoImage+=0.06
albedoImage[w]=-32768.



#load the shape model

###use the following as reference for plotting
#http://stackoverflow.com/questions/18897786/transparency-for-poly3dcollection-plot-in-matplotlib
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

#get the vertices indices
vertIndices=[]
for i in range(6*n**2+2,len(data)):
    s=data[i].split()
    vertIndices.append([int(float(s[1]))-1,int(float(s[2]))-1,int(float(s[3]))-1])
vertIndices=num.array(vertIndices)

vertlons=(num.arctan2(vertices[:,1],vertices[:,0])*r2d)#%360
vertlats=num.arcsin(vertices[:,2]/(vertices[:,0]**2+vertices[:,1]**2+vertices[:,2]**2)**0.5)*r2d
vertLl=num.zeros((len(vertlons),2)).astype('float64')
vertLl[:,0]=vertlons
vertLl[:,1]=vertlats
vivll=vertLl[vertIndices]

whichVertsToPlot=[]
for i in range(len(vivll)):
    if num.max(vivll[i,:,0])-num.min(vivll[i,:,0])<300:
        whichVertsToPlot.append(i)
whichVertsToPlot=num.array(whichVertsToPlot)
#whichVertsToPlot=num.arange(len(vivll))

"""
#testing
cmap=pyl.get_cmap('jet')
lLCollection=PolyCollection(vivll[whichVertsToPlot])
collectionColours=cmap(num.ones(len(vivll)))
collectionColours[:,3]=0.8
lLCollection.set_facecolors(collectionColours[whichVertsToPlot])
lLCollection.set_linewidths(0.0)


fig=pyl.figure()
sp=fig.add_subplot(111)
ax=sp.add_collection(lLCollection)
cheat=sp.scatter([-10,-10],[100,100],c=num.array([0,1]))
pyl.colorbar(cheat)
sp.set_xlim(0,360)
sp.set_ylim(-90,90)
pyl.show()
sys.exit()
#end testing
"""


#get the vertices indices
(mids,normals)=midsNormals(vertices[vertIndices])
lons=(num.arctan2(mids[:,1],mids[:,0])*r2d)#%360
lats=num.arcsin(mids[:,2]/(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5)*r2d
D=(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5
meanRadius=num.median(D)
print 'Median radius is:',meanRadius


gridLong=num.linspace(num.min(lons),num.max(lons),360*10)
gridLat=num.linspace(num.min(lats),num.max(lats),360*10)
gridD=interp.griddata((lons,lats),D,(gridLong[None,:],gridLat[:,None]),method='linear')

#this is to put lons in the same zero as plotShape_VIMS_doubleArm.py
lonsNotMod=lons*1.
lons=lons%360

albedos=getAlbedo(lonsNotMod,lats,albedoImage)


resolutions=[]
imageNames=['2004163T121836_2004163T192848/cv1465649433_1',
'2004163T121836_2004163T192848/cv1465649746_1',
'2004163T121836_2004163T192848/cv1465649834_1',
'2004163T121836_2004163T192848/cv1465649979_1',
'2004163T121836_2004163T192848/cv1465650070_1',
'2004163T121836_2004163T192848/cv1465650745_1',
'2004163T121836_2004163T192848/cv1465650834_1',
'2004163T121836_2004163T192848/cv1465651001_1',
'2004163T121836_2004163T192848/cv1465651857_1',
'2004163T121836_2004163T192848/cv1465661929_1',
'2004163T121836_2004163T192848/cv1465662167_1',
'2004163T121836_2004163T192848/cv1465662631_1',
'2004163T121836_2004163T192848/cv1465664774_1',
'2004163T121836_2004163T192848/cv1465669068_1',
'2004163T121836_2004163T192848/cv1465669944_1',
'2004163T121836_2004163T192848/cv1465670212_1',
'2004163T121836_2004163T192848/cv1465670650_1',
'2004163T121836_2004163T192848/cv1465672161_1',
'2004163T121836_2004163T192848/cv1465673600_1',
'2004163T193015_2004164T051726/cv1465677443_1',
'2004163T193015_2004164T051726/cv1465677670_1',
'2004163T193015_2004164T051726/cv1465678419_1',
'2004163T193015_2004164T051726/cv1465678911_1',
'2004163T193015_2004164T051726/cv1465679413_1',
'2004163T193015_2004164T051726/cv1465679675_1',
'2004163T193015_2004164T051726/cv1465679932_1']


initLimit=5.

extractSpec=False

for imageName in imageNames:
    if not extractSpec:continue
    waterDepths=[]
    siliDepths=[]
    for jj in range(len(mids)):
        waterDepths.append([])
        siliDepths.append([])
        
    print 'Working with image:'
    print imageName


    (visDataMin,IRDataMin,oxv,oyv,oxi,oyi,distancenot)=getPars('/data/VIMS/covims_0004/procdata/%s.parFile'%(imageName))
    offsetsVis=num.array([oxv,oyv])
    offsetsIR=num.array([oxi,oyi])

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
    IRBGPixels=num.where(imDataIR<256.)
    imDataIR[IRBGPixels]=0.0

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

    imData=num.concatenate([[imDataVis],[imDataIR],[maskVis]])
    (dump,A,B)=imData.shape
    
    
    with open('/data/VIMS/covims_0004/procdata/%s_vis.campt.pickle'%(imageName)) as han:
        latLongObjVis=pick.load(han)
    
    with open('/data/VIMS/covims_0004/procdata/%s_ir.campt.pickle'%(imageName)) as han:
        latLongObjIR=pick.load(han)
    
    
    
    
    
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
    resolutions.append(sampResolutionsIR[argI])
    
    
    
    
    w=[[0]]
    inot=int(CassIR[w][:,5])#same as mid_i
    jnot=int(CassIR[w][:,6])#same as mid_j
    Tnot=CassIR[w][0][3]
    
    
    
    pixelTimesVis=getPixTimes(CassVis,imDataVis.shape[0],imDataVis.shape[1],Tnot,inot,jnot)
    pixelTimesIR=getPixTimes(CassIR,imDataIR.shape[0],imDataIR.shape[1],Tnot,inot,jnot)


    #distance units are in km, resolution units in km/pix
    [Xnot,Ynot,Znot]=CassIR[w][0][:3]
    if distancenot==None:
        distancenot=CassIR[w][0][4]
    else:
        Xnot*=distancenot/CassIR[w][0][0]
        Ynot*=distancenot/CassIR[w][0][0]
        Znot*=distancenot/CassIR[w][0][0]
    
    
    #nominal values from the VIMS manual are 0.17*3 and 0.5 mrad
    resRat=num.min(sampResolutionsVis[W])/sampResolutionsIR[argI]
    if resRat<0.5:
        pixScaleVis=0.17
    else:
        pixScaleVis=0.17*3
    pixScaleIR=0.5
    
    (junk,deltaXVis,deltaYVis,deltaZVis)=getVels(CassVis,pixelTimesVis,Xnot,Ynot,Znot,inot,jnot,polyOrder=2)
    (spaceCraftVectornot,deltaXIR,deltaYIR,deltaZIR)=getVels(CassIR,pixelTimesIR,Xnot,Ynot,Znot,inot,jnot,polyOrder=2)
    
    long_o_not=latLongObjIR[inot][jnot]['SubSpacecraftLongitude']
    lat_o_not=latLongObjIR[inot][jnot]['SubSpacecraftLatitude']
    az_o_not=latLongObjIR[inot][jnot]['SpacecraftAzimuth']
    long_s=latLongObjIR[inot][jnot]['SubSolarLongitude']
    lat_s=latLongObjIR[inot][jnot]['SubSolarLatitude']
    az_s=latLongObjIR[inot][jnot]['SubSolarAzimuth']
    initial=num.array([long_o_not,lat_o_not,az_o_not,long_s,lat_s])


    #load the best point
    (bestPoint,goodSamps)=getFit('/data/VIMS/covims_0004/procdata/%s.fit_pickle'%(imageName),limit=initLimit,initialAngles=initial)
    [long_o_not,lat_o_not,az_o_not,long_s,lat_s,distancenot,offXV,offYV,offXI,offYI,offsetVel,offsetVelAngle,chi]=bestPoint
    print 'EMCEE best chi:',chi

    gotten=False
    if os.path.isfile('/data/VIMS/covims_0004/procdata/%s.rand'%(imageName)):
        with open('/data/VIMS/covims_0004/procdata/%s.rand'%(imageName)) as randHan:
            randoms=randHan.readlines()
        for i in range(0,len(randoms),4):
            s=''
            for j in range(4):
                s+=randoms[i+j].split('\n')[0]
            s=s.split()
            rchi=float(s[len(s)-1])
            if rchi>chi or not gotten:
                x=[]
                for j in range(1,len(s)-1):
                    x.append(float(s[j].split(']')[0]))
                x=num.array(x)
                [long_o_not,lat_o_not,az_o_not,long_s,lat_s,distancenot,offXV,offYV,offXI,offYI,offsetVel,offsetVelAngle]=x
                chi=rchi
                gotten=True
        print 'RANDOM best chi:',chi
    offsetsVis=num.array([offXV,offYV])
    offsetsIR =num.array([offXI,offYI])

    sampleResolutionnot=pixScaleIR*distancenot/1000.
    print ' IR fitted nominal sample resolution:',sampleResolutionnot
    print ' IR header sample resolution:',sampResolutionsIR[argI]/1000.

    #print 'Vis nominal sample resolution:',num.min(sampResolutionsVis[W])
    
    #calculate the shape model
    #NEED az_adjust=0 to get the pixels in the right places!
    (imageVis,poly3d,colours,rot_vert,vertsInPixelVis,chiv)=shapeGen_VIMS(vertices,vertIndices,
                                                    long_o_not,lat_o_not,az_o_not,
                                                    long_s,lat_s,0,
                                                    distancenot,
                                                    offsetsVis, offsetVel,offsetVelAngle,
                                                    pixScaleVis,
                                                    inot,jnot,
                                                    deltaXVis,deltaYVis,deltaZVis,
                                                    lons,lats,
                                                    pixelTimesVis,
                                                    imData[0],vis=True,mask=imData[2],az_adjust=-90)
    (imageIR,poly3d,colours,rot_vert,vertsInPixelIR,chii)=shapeGen_VIMS(vertices,vertIndices,
                                                    long_o_not,lat_o_not,az_o_not,
                                                    long_s,lat_s,0,
                                                    distancenot,
                                                    offsetsIR, offsetVel,offsetVelAngle,
                                                    pixScaleIR,
                                                    inot,jnot,
                                                    deltaXIR,deltaYIR,deltaZIR,
                                                    lons,lats,
                                                    pixelTimesVis,
                                                    imData[1],vis=False,az_adjust=-90)


    #produce the spectra files
    with pyf.open('/data/VIMS/covims_0004/procdata/'+imageName+'_ir.fits') as shan:
        IRspecData=shan[0].data
    with pyf.open('/data/VIMS/covims_0004/procdata/'+imageName+'_vis_strpd.fits') as shan:
        VISspecData=shan[0].data

    if '2004163T193015_2004164T051726' in imageName:
        print 'no spectral swap'
    else:
        IRspecData=IRspecData[:,::-1,::-1]
        VISspecData=VISspecData[:,::-1,::-1]
    (IRchannels,A,B)=IRspecData.shape
    (VISchannels,a,b)=VISspecData.shape


    IRSpec=[]
    IRvip=[]
    for l in range(A):
        for s in range(B):
            (w,spec)=specAnalysis.getSpec(IRspecData,l=l,s=s)
            med=num.nanmedian(spec)
            if len(vertsInPixelIR[l][s])>0:#med>IRDataMin:# and imData[1][l][s]>0:
                IRSpec.append(spec)
                IRvip.append(vertsInPixelIR[l][s])

    VISSpec=[]
    VISvip=[]
    for l in range(a):
        for s in range(b):
            (w,spec)=specAnalysis.getSpecVis(VISspecData,l=l,s=s)
            med=num.nanmedian(spec)
            #print med,visDataMin,len(vertsInPixelIR[l][s]),imData[0][l][s]
            if len(vertsInPixelVis[l][s])>0:#med>visDataMin:# and imData[0][l][s]>0:
                VISSpec.append(spec)
                VISvip.append(vertsInPixelVis[l][s])
    if extractSpec:
        with open('/data/VIMS/covims_0004/procdata/'+imageName+'_spec.pickle','w+') as han:
            pick.dump([VISSpec,IRSpec,VISvip,IRvip,sampleResolutionnot],han)
            #sys.exit()

################
####plot parameters
cmapToUse='YlOrRd' #jet

generateMeanSpectra=False
dumpingRegions=False
saveMap=False
showBinaryColours=False
showContours=True
showMap=True
showJason=False
showScatter=False
showShape=False
saveMovie=False
includeMedLow=True
threeColour=False
doLineProfile=True

cmap=pyl.get_cmap(cmapToUse)

numContours=25
if includeMedLow: minPix=1
else: minPix=1

if showMap:
    if showJason:
        fs=(10,10)
        alpha=1.
    else:
        fs=(10,10)
        alpha=1.
else:
    fs=(15,10)
    alpha=1.


whichRes='low'

if whichRes=='high':
    #high res
    resLow=0
    resHigh=20.
    title='High Resolution, <16 km/pix'
elif whichRes=='med':
    #med res
    resLow=20.
    resHigh=50.
    title='Medium Resolution, 22>r>16 km/pix'
elif whichRes=='low':
    #low res
    resLow=50.
    resHigh=100
    title='Low Resolution, >22 km/pix'


#0.1 to 0.41 work well for the global map
cMaxw=0.41
cMinw=0.1
cbTitlew='1.6+2 Water-ice Absorption'

cMaxwb=0.81
cMinwb=0.64
cbTitlewb='3 Water-ice Absorption'

cMaxww=0.5
cMinww=0.07
cbTitleww='Absorption Ratio'


cMaxo=-0.03
cMino=-0.09
cbTitleo='R-band Reddening'

cMind=-0.0
cMaxd=0.16
cbTitled='B-band Reddening'
####
#################
extent=[-180,180,-90,90]


if dumpingRegions:
    dumpingSpectra=[[],[],[],[]]
if generateMeanSpectra:
    whichSpecVis=[]
    whichSpecIR=[]
    for i in range(len(mids)):
        whichSpecVis.append([])
        whichSpecIR.append([])


    VISSpecs=[]
    IRSpecs=[]
    for imageName in imageNames:
        with open('/data/VIMS/covims_0004/procdata/'+imageName+'_spec.pickle') as han:
            [VISSpec,IRSpec,VISvip,IRvip,sampleResolutionnot]=pick.load(han)
            if sampleResolutionnot<resLow or sampleResolutionnot>resHigh: continue
            print imageName, ' with resolution',sampleResolutionnot
            for i in range(len(VISSpec)):
                if len(VISvip[i])==0:continue
                VISSpecs.append(VISSpec[i]*1.0)
                l=len(VISSpecs)-1
                for j in VISvip[i]:
                    whichSpecVis[int(j)].append(l)
            for i in range(len(IRSpec)):
                if len(IRvip[i])==0:continue
                IRSpecs.append(IRSpec[i]*1.0)
                l=len(IRSpecs)-1
                for j in IRvip[i]:
                    whichSpecIR[int(j)].append(l)


    plotted=0
    avSpecs=[]
    for i in range(len(whichSpecIR)):
        v=[]
        ir=[]
        noVisSpec=True
        if len(whichSpecIR[i])>=1:
            for l in whichSpecVis[i]:
                s=normVisSpec(VISSpecs[l])
                v.append(s)
            if v<>[]:
                v=num.array(v)
                noVisSpec=False
            else: v=[1.0]
            x=num.median(v)
            v/=x
            vs=[]
            for l in whichSpecVis[i]:
                vs.append(VISSpecs[l]/v[len(vs)])
            if vs<>[]:
                vs=num.array(vs)
            else:
                vs=num.ones([1,len(viswave)]).astype('float')


            for l in whichSpecIR[i]:
                s=normIRSpec(IRSpecs[l])
                ir.append(s)
            ir=num.array(ir)
            x=num.median(ir)
            ir/=x
            irs=[]
            for l in whichSpecIR[i]:
                #pyl.plot(irwave,IRSpecs[l]/ir[len(irs)])
                irs.append(IRSpecs[l]/ir[len(irs)])
            #print 'Do double check that the way you are medianing these spectra is intelligent.'

            medVis=nanmean(num.array(vs))
            medVis/=normVisSpec(medVis)
            junk=num.copy(medVis)
            medIR=nanmean(num.array(irs))
            medIR/=normIRSpec(medIR)

            w=num.where(irwave<viswave[86])
            f=interp.interp1d(irwave[w],medIR[w])
            w=num.where(viswave[:86]>irwave[0])
            medVis[w]=(medVis[w]+f(viswave[w]))/2.
            w=num.where(irwave>=viswave[86])
            fullspec=num.concatenate([medVis[:86],medIR[w]])
            fullwave=num.concatenate([viswave[:86],irwave[w]])


            #if len(irs)>=2:
            #    pyl.plot(fullwave,fullspec)
            #    pyl.plot(irwave,IRSpecs[0])
            #    pyl.plot(irwave,IRSpecs[1])
            #    pyl.show()
            #    sys.exit()
            #take the median of water depths, not the water depth from the median spectrum
            nSpec=len(irs)
            wats=[]
            watbs=[]
            for l in range(nSpec):
                (wat,watb)=specAnalysis.water(irwave,irs[l])
                wats.append(wat)
                watbs.append(watb)

            wat,watb=num.median(wats),num.median(watbs)

            if not noVisSpec:
                oSlope=specAnalysis.oSlope(fullwave,fullspec)
                dSlope=specAnalysis.dSlope(fullwave,fullspec)
            else:
                oSlope,dSlope=-32768.,-32768.

            if dumpingRegions and whichRes=='high':
                if abs(wat-0.1)<0.03 and abs(watb-0.66)<0.015: #low
                    dumpingSpectra[0].append(num.copy(fullspec))
                elif abs(wat-0.21)<0.03 and abs(watb-0.73)<0.03: #mid
                    dumpingSpectra[1].append(num.copy(fullspec))
                elif abs(wat-0.30)<0.03 and abs(watb-0.78)<0.03: #high
                    dumpingSpectra[2].append(num.copy(fullspec))
                elif wat>0.37 and watb>0.77 and watb<0.82: #rich
                    dumpingSpectra[3].append(num.copy(fullspec))

            #print wat
            #if abs(wat-.29)<0.05 and abs(watb-0.78)<0.03:
                #pyl.clf()
                #pyl.plot(fullwave,fullspec,'k-')
                #pyl.title(str(wat)+' '+str(watb))
                #dumpingSpectra.append(num.copy(fullspec))
                #pyl.show()
                #sys.exit()

            plotted+=1
            #i==44331 in the 'med' is a problem spectrum
            if plotted in range(0,22000,1000):
                pyl.plot(viswave[:86],junk[:86],lw=2)
                pyl.plot(irwave,medIR,lw=2)
                pyl.plot(fullwave,fullspec)

            avSpecs.append([fullspec,len(whichSpecVis[i]),len(whichSpecIR[i]),wat,watb,oSlope,dSlope])


        else:
            avSpecs.append([None,len(whichSpecVis[i]),len(whichSpecIR[i]),None,None,None,None])

    #print 'number of specs',len(avSpecs)
    #pyl.xlabel('$\\lambda \\mbox{ ($\\mu$m)}$')
    #pyl.ylabel('Normalized Reflectance')
    #pyl.show()
    #sys.exit()
    #print plotted

    if dumpingRegions and whichRes=='high':
        with open('WaterModelling/lowSpec_highres.pickle','w+') as outhan:
            pick.dump([fullwave,dumpingSpectra[0]],outhan)
        with open('WaterModelling/midSpec_highres.pickle','w+') as outhan:
            pick.dump([fullwave,dumpingSpectra[1]],outhan)
        with open('WaterModelling/highSpec_highres.pickle','w+') as outhan:
            pick.dump([fullwave,dumpingSpectra[2]],outhan)
        with open('WaterModelling/richSpec_highres.pickle','w+') as outhan:
            pick.dump([fullwave,dumpingSpectra[3]],outhan)
    with open(whichRes+'_avSpecs.pickle','w+') as han:
        pick.dump(avSpecs,han)
    pyl.show()
    #exit()
    sys.exit()
avSpecs={}
with open('high_avSpecs.pickle') as han:
    avSpecs['high']=pick.load(han)
with open('med_avSpecs.pickle') as han:
    avSpecs['med']=pick.load(han)
with open('low_avSpecs.pickle') as han:
    avSpecs['low']=pick.load(han)

specs={'high':[],'med':[],'low':[]}
numVisSamps={'high':[],'med':[],'low':[]}
numIRSamps={'high':[],'med':[],'low':[]}
waterDepths={'high':[],'med':[],'low':[]}
boundWaterDepths={'high':[],'med':[],'low':[]}
oSlopes={'high':[],'med':[],'low':[]}
dSlopes={'high':[],'med':[],'low':[]}
for k in avSpecs:
    for i in range(len(avSpecs[k])):
        specs[k].append(avSpecs[k][i][0])
        numVisSamps[k].append(avSpecs[k][i][1])
        numIRSamps[k].append(avSpecs[k][i][2])
        #avSpecs.append([fullspec,len(whichSpecVis[i]),len(whichSpecIR[i]),wat,watb,oSlope,dSlope])
        if avSpecs[k][i][3]==None or avSpecs[k][i][2]<minPix or not isValidSpecValues(avSpecs[k][i]):
            waterDepths[k].append(-32768.)
            boundWaterDepths[k].append(-32768.)
            oSlopes[k].append(-32768.)
            dSlopes[k].append(-32768.)
        else:
            waterDepths[k].append(avSpecs[k][i][3])
            boundWaterDepths[k].append(avSpecs[k][i][4])
            oSlopes[k].append(avSpecs[k][i][5])
            dSlopes[k].append(avSpecs[k][i][6])
    specs[k]=num.array(specs[k])
    numVisSamps[k]=num.array(numVisSamps[k])
    numIRSamps[k]=num.array(numIRSamps[k])
    waterDepths[k]=num.array(waterDepths[k])
    boundWaterDepths[k]=num.array(boundWaterDepths[k])
    oSlopes[k]=num.array(oSlopes[k])
    dSlopes[k]=num.array(dSlopes[k])


if showScatter:
    elevation=(mids[:,0]**2+mids[:,1]**2+mids[:,2]**2)**0.5
    w=num.where((waterDepths['high']<>-32768.)&(boundWaterDepths['high']<>-32768.)&(oSlopes['high']<>-32768.))#&(oSlopes['high']<0.02)&(dSlopes['high']<>-32768.))


    ww=num.where(albedos[w]>=0.06)
    pyl.scatter(albedos[w][ww],waterDepths['high'][w][ww],alpha=0.1,c='k')
    pyl.xlabel('Visual Albedo')
    pyl.ylabel('1.5+2~$\mu$m Absorption Depth')
    pyl.show()
    #W=num.where((waterDepths['high'][w]<0.34)&(boundWaterDepths['high'][w]>0.75))
    #pyl.scatter(lonsNotMod[w][W],lats[w][W])
    #pyl.show()



    pyl.scatter(waterDepths['high'][w],boundWaterDepths['high'][w],c='k',alpha=0.1)
    x=num.array([0.05,0.4])
    y=0.6+0.6*x
    pyl.plot(x,y,'r-',lw=2)
    delt=boundWaterDepths['high'][w]-(0.6*waterDepths['high'][w]+0.6)
    pyl.scatter(waterDepths['high'][w],delt,c='r',alpha=0.2)
    pyl.xlabel('$1.55+2 \\mbox{ $\\mu$m}$ absorption')
    pyl.ylabel('$3 \\mbox{ $\\mu$m}$ absorption')
    pyl.show()
    sys.exit()
    #fig3d=pyl.figure('3dScatter')
    #ax=fig3d.add_subplot(111,projection='3d')
    #ax.scatter(waterDepths['high'][w],oSlopes['high'][w],dSlopes['high'][w],marker='.',alpha=0.6)
    #ax.set_xlabel('Water-Ice Absorption')
    #ax.set_ylabel('R-band Slope')
    #ax.set_zlabel('B-Band Slope')

    specFig=pyl.figure('JasonSpec')
    
    woUniquesInJason=[]
    elevUniquesInJason=[]
    for i in range(len(waterDepths['high'][w])):
        x=[waterDepths['high'][w][i],oSlopes['high'][w][i],dSlopes['high'][w][i]]
        if x not in woUniquesInJason:
            W=num.where((waterDepths['high'][w]==x[0])&(oSlopes['high'][w]==x[1])&(dSlopes['high'][w]==x[2]))
            meanLon=num.mean(lonsNotMod[w][W])
            meanLat=num.mean(lats[w][W])
            meanElev=num.mean(elevation[w][W])
            if (meanLon>-5 and meanLon<65 and meanLat>-15 and meanLat<50) and not (meanLon>19 and meanLon<30 and meanLat>9 and meanLat<21):
                woUniquesInJason.append(x)
                elevUniquesInJason.append(meanElev)

                if meanElev>106.:
                    pyl.plot(specs['high'][w][i],'r-')
    #pyl.show()
    #sys.exit()

    woUniquesOutJason=[]
    elevUniquesOutJason=[]
    for i in range(len(waterDepths['high'][w])):
        x=[waterDepths['high'][w][i],oSlopes['high'][w][i],dSlopes['high'][w][i]]
        if x not in woUniquesOutJason:
            W=num.where((waterDepths['high'][w]==x[0])&(oSlopes['high'][w]==x[1])&(dSlopes['high'][w]==x[2]))
            meanLon=num.mean(lonsNotMod[w][W])
            meanLat=num.mean(lats[w][W])
            meanElev=num.mean(elevation[w][W])
            if (meanLon>-60 and meanLon<-5 and meanLat>-15 and meanLat<50):
    
                woUniquesOutJason.append(x)
                elevUniquesOutJason.append(meanElev)


    
    woUniquesFloorJason=[]
    elevUniquesFloorJason=[]
    for i in range(len(waterDepths['high'][w])):
        x=[waterDepths['high'][w][i],oSlopes['high'][w][i],dSlopes['high'][w][i]]
        if x not in woUniquesFloorJason:
            W=num.where((waterDepths['high'][w]==x[0])&(oSlopes['high'][w]==x[1])&(dSlopes['high'][w]==x[2]))
            meanLon=num.mean(lonsNotMod[w][W])
            meanLat=num.mean(lats[w][W])
            meanElev=num.mean(elevation[w][W])
            if (meanLon>19 and meanLon<30 and meanLat>9 and meanLat<21):
                woUniquesFloorJason.append(x)
                elevUniquesFloorJason.append(meanElev)

                pyl.plot(specs['high'][w][i],'k-')
    pyl.show()
    sys.exit()
    
    woUniquesInJason=num.array(woUniquesInJason)
    woUniquesOutJason=num.array(woUniquesOutJason)
    woUniquesFloorJason=num.array(woUniquesFloorJason)

    fige=pyl.figure('elev',figsize=(10,8))
    fige.subplots_adjust(hspace=0)
    sp1=fige.add_subplot(311,xticklabels='')
    pyl.scatter(elevUniquesInJason,woUniquesInJason[:,0],c='b',alpha=0.3)
    pyl.scatter(elevUniquesOutJason,woUniquesOutJason[:,0],c='r',alpha=0.3)
    pyl.scatter(elevUniquesFloorJason,woUniquesFloorJason[:,0],c='y',alpha=0.7)
    pyl.ylabel('Water-ice Absorption')
    print 'Elevation water-ice inside Jason',stats.spearmanr(elevUniquesInJason,woUniquesInJason[:,0])
    print 'Elevation water-ice out of Jason',stats.spearmanr(elevUniquesOutJason,woUniquesOutJason[:,0])
    print 'Elevation water-ice floor  Jason',stats.spearmanr(elevUniquesFloorJason,woUniquesFloorJason[:,0])
    print
    sp2=fige.add_subplot(312)
    pyl.scatter(elevUniquesInJason,woUniquesInJason[:,1],c='b',alpha=0.3)
    pyl.scatter(elevUniquesOutJason,woUniquesOutJason[:,1],c='r',alpha=0.3)
    pyl.scatter(elevUniquesFloorJason,woUniquesFloorJason[:,1],c='y',alpha=0.7)
    pyl.xlabel('Elevation (km)')
    pyl.ylabel('R-band Slope (\% per 100 nm)')
    print 'Elevation oSlope inside Jason',stats.spearmanr(elevUniquesInJason,woUniquesInJason[:,1])
    print 'Elevation oSlope out of Jason',stats.spearmanr(elevUniquesOutJason,woUniquesOutJason[:,1])
    print 'Elevation oSlope floor  Jason',stats.spearmanr(elevUniquesFloorJason,woUniquesFloorJason[:,1])
    print
    sp3=fige.add_subplot(313)
    pyl.scatter(elevUniquesInJason,woUniquesInJason[:,2],c='b',alpha=0.3)
    pyl.scatter(elevUniquesOutJason,woUniquesOutJason[:,2],c='r',alpha=0.3)
    pyl.scatter(elevUniquesFloorJason,woUniquesFloorJason[:,2],c='y',alpha=0.7)
    pyl.xlabel('Elevation (km)')
    pyl.ylabel('B-band Slope (\% per 100 nm)')
    print 'Elevation dSlope inside Jason',stats.spearmanr(elevUniquesInJason,woUniquesInJason[:,2])
    print 'Elevation dSlope out of Jason',stats.spearmanr(elevUniquesOutJason,woUniquesOutJason[:,2])
    print 'Elevation dSlope floor  Jason',stats.spearmanr(elevUniquesFloorJason,woUniquesFloorJason[:,2])
    print
    
    figCorr=pyl.figure('Correlations')
    print 'Water-ice oSlope inside Jason',stats.spearmanr(woUniquesInJason[:,0],woUniquesInJason[:,1])
    print 'Water-ice oSlope out of Jason',stats.spearmanr(woUniquesOutJason[:,0],woUniquesOutJason[:,1])
    print 'Water-ice oSlope Floor Jason',stats.spearmanr(woUniquesFloorJason[:,0],woUniquesFloorJason[:,1])
    print
    print 'oSlope dSlope inside Jason',stats.spearmanr(woUniquesInJason[:,2],woUniquesInJason[:,1])
    print 'oSlope dSlope out of Jason',stats.spearmanr(woUniquesOutJason[:,2],woUniquesOutJason[:,1])
    print 'oSlope dSlope Floor Jason',stats.spearmanr(woUniquesFloorJason[:,2],woUniquesFloorJason[:,1])
    pyl.scatter(woUniquesInJason[:,0],woUniquesInJason[:,1],marker='o',alpha=0.3,c='b')
    pyl.scatter(woUniquesOutJason[:,0],woUniquesOutJason[:,1],marker='o',alpha=0.3,c='r')
    pyl.scatter(woUniquesFloorJason[:,0],woUniquesFloorJason[:,1],marker='o',alpha=1.,c='y')
    pyl.xlabel('Water Absorption Depth')
    pyl.ylabel('0.5-0.7 micron Optical Slope (\% per 100 nm)')
    pyl.show()
    sys.exit()



fullWaterDepths=waterDepths['high']*1.0
if includeMedLow:
    w=num.where(fullWaterDepths==-32768.)
    for i in w[0]:
        fullWaterDepths[i]=waterDepths['med'][i]
    w=num.where(fullWaterDepths==-32768.)
    for i in w[0]:
        fullWaterDepths[i]=waterDepths['low'][i]

fullBoundWaterDepths=boundWaterDepths['high']*1.0
if includeMedLow:
    w=num.where(fullBoundWaterDepths==-32768.)
    for i in w[0]:
        fullBoundWaterDepths[i]=boundWaterDepths['med'][i]
    w=num.where(fullBoundWaterDepths==-32768.)
    for i in w[0]:
        fullBoundWaterDepths[i]=boundWaterDepths['low'][i]


fulloSlopes=oSlopes['high']*1.0
if includeMedLow:
    w=num.where((fulloSlopes==-32768.)|(fulloSlopes>0.02))
    for i in w[0]:
        fulloSlopes[i]=oSlopes['med'][i]
    w=num.where((fulloSlopes==-32768.)|(fulloSlopes>0.02))
    for i in w[0]:
        fulloSlopes[i]=oSlopes['low'][i]

fulldSlopes=dSlopes['high']*1.0
if includeMedLow:
    w=num.where((fulldSlopes==-32768.))
    for i in w[0]:
        fulldSlopes[i]=dSlopes['med'][i]
    w=num.where((fulldSlopes==-32768.))
    for i in w[0]:
        fulldSlopes[i]=dSlopes['low'][i]


colw=num.zeros(len(fullWaterDepths)).astype('float')
alphasw=colw*0.0
whichVertsToColour=num.where(fullWaterDepths<>-32768.)
colw[whichVertsToColour]=(fullWaterDepths[whichVertsToColour]-cMinw)/(cMaxw-cMinw)
alphasw[whichVertsToColour]=alpha


lLCollection1=PolyCollection(vivll[whichVertsToPlot],zorder=10)
collectionColoursw=cmap(colw)
collectionColoursw[:,3]=alphasw
#collectionColours[w][:,3]=0.0
lLCollection1.set_facecolors(collectionColoursw[whichVertsToPlot])
lLCollection1.set_linewidths(0.0)

colwb=num.zeros(len(fullBoundWaterDepths)).astype('float')
alphaswb=colwb*0.0
whichVertsToColour=num.where(fullBoundWaterDepths<>-32768.)
colwb[whichVertsToColour]=(fullBoundWaterDepths[whichVertsToColour]-cMinwb)/(cMaxwb-cMinwb)
alphaswb[whichVertsToColour]=alpha

lLCollection4=PolyCollection(vivll[whichVertsToPlot],zorder=10)
collectionColourswb=cmap(colwb)
collectionColourswb[:,3]=alphaswb
#collectionColours[w][:,3]=0.0
lLCollection4.set_facecolors(collectionColourswb[whichVertsToPlot])
lLCollection4.set_linewidths(0.0)


colww=num.zeros(len(fullBoundWaterDepths)).astype('float')
alphasww=colww*0.0
wwRat=fullWaterDepths/fullBoundWaterDepths
#wwRat[num.where(num.isnan(wwRat))]=0.0
#wwRat[num.where(( fullWaterDepths<0.1) | (fullBoundWaterDepths<0.6))]=0.0
#whichVertsToColour=num.where(( fullWaterDepths>0.01)& (fullBoundWaterDepths>0.6))
#colww[whichVertsToColour]=(wwRat[whichVertsToColour]-cMinww)/(cMaxww-cMinww)
w=num.where((fullWaterDepths>0.36) & (fullBoundWaterDepths>0.75))
colww[w]=1.
w=num.where((fullWaterDepths<0.16))
colww[w]=0.1
w=num.where((fullWaterDepths>=0.16)&(fullWaterDepths<0.26))
colww[w]=0.4
w=num.where((fullWaterDepths>=0.26)&(fullWaterDepths<0.36))
colww[w]=0.7
alphasww[whichVertsToColour]=alpha

lLCollection5=PolyCollection(vivll[whichVertsToPlot],zorder=10)
collectionColoursww=cmap(colww)
collectionColoursww[:,3]=alphasww
#collectionColours[w][:,3]=0.0
lLCollection5.set_facecolors(collectionColoursww[whichVertsToPlot])
lLCollection5.set_linewidths(0.0)
lLCollection5.set_linewidths(0.0)


colo=num.zeros(len(fullWaterDepths)).astype('float')
alphaso=colw*0.0
whichVertsToColour=num.where((fulloSlopes>-32768.)&(fulloSlopes<0.2))
colo[whichVertsToColour]=(fulloSlopes[whichVertsToColour]-cMino)/(cMaxo-cMino)
alphaso[whichVertsToColour]=alpha


lLCollection2=PolyCollection(vivll[whichVertsToPlot],zorder=10)
collectionColourso=cmap(colo)
collectionColourso[:,3]=alphaso
#collectionColours[w][:,3]=0.0
lLCollection2.set_facecolors(collectionColourso[whichVertsToPlot])
lLCollection2.set_linewidths(0.0)


cold=num.zeros(len(fullWaterDepths)).astype('float')
alphasd=colw*0.0
whichVertsToColour=num.where((fulldSlopes>-32768.))
cold[whichVertsToColour]=(fulldSlopes[whichVertsToColour]-cMind)/(cMaxd-cMind)
alphasd[whichVertsToColour]=alpha

lLCollection3=PolyCollection(vivll[whichVertsToPlot],zorder=10)
collectionColoursd=cmap(cold)
collectionColoursd[:,3]=alphasd
#collectionColours[w][:,3]=0.0
lLCollection3.set_facecolors(collectionColoursd[whichVertsToPlot])
lLCollection3.set_linewidths(0.0)

if threeColour:

    mixColours=num.ones((len(colo),4)).astype('float64')
    mixColours[:,0]=colw
    mixColours[:,1]=colwb
    mixColours[:,2]=0.5
    mixColours[:,3]=1.0
    mixColours=num.clip(mixColours,0.,1.)
    w=num.where((mixColours[:,0]==0)&(mixColours[:,1]==0)&(mixColours[:,2]==0))
    mixColours[:,3][w]=0.0

    mixColours=cmap(colww)

    
    lLCollection3C=PolyCollection(vivll[whichVertsToPlot],zorder=10)
    lLCollection3C.set_facecolors(mixColours[whichVertsToPlot])
    lLCollection3C.set_linewidths(0.0)

    fig3c=pyl.figure('threeColour',figsize=(10,8))
    fig3c.subplots_adjust(hspace=0)
    
    sp1=fig3c.add_subplot(111)
    #sp1.imshow(newImage,interpolation='nearest',cmap='gray',extent=extent,zorder=1)
    sp1.add_collection(lLCollection3C)
    sp1.set_xlim(-80,120)
    sp1.set_ylim(-40,60)


    if showContours:
        contours=pyl.contour(gridLong,gridLat,gridD-106,numContours,colors='w',zorder=10)
        contours.levels = [nf(val) for val in contours.levels ]
        # Label levels with specially formatted floats
        if pyl.rcParams["text.usetex"]:
            fmt = r'%r'
        else:
            fmt = '%r'
        pyl.clabel(contours, contours.levels, inline=True, fmt=fmt, fontsize=14)


    if doLineProfile:
        lp=profileHandler(fig3c,sp1,lonsNotMod,lats,D,mids,mixColours,mixColours=mixColours)
        pyl.show()

    if doLineProfile: sys.exit()
    
if showShape:
    #if not saveMovie:
    shapeFig=pyl.figure('Shape',figsize=(15,8))
    shapeFig.subplots_adjust(hspace=0,wspace=0)
    lonlon=num.array([0.,72.,144.,216.,288.])+35.
    latlat=[35.,-35.]
    for i in range(len(lonlon)):
        for k in range(len(latlat)):
            ax1=shapeFig.add_subplot(len(latlat),len(lonlon),
                                     1+k*len(lonlon)+i,
                                     projection='3d')
            ax1.set_aspect('equal')
            ax1.set_xlim(-120,120)
            ax1.set_ylim(-120,120)
            ax1.set_zlim(-120,120)
            ax1.set_axis_bgcolor('0.0')
            ax1.view_init(azim=0, elev=0)
            ax1.set_axis_off()
            pyl.title('H$_2$O Absorption',color='w')

            poly3d=reshape(vertices,vertIndices,
            lonlon[i],latlat[k],0.,
            80.,-12.,0.)[1]

            w=num.where((collectionColoursw[:,0]==1.)&(collectionColoursw[:,1]==1.))[0]
            for j in w:
                collectionColoursw[j]=num.array([0.05,0.05,0.05,1.0])
            w=num.where((collectionColourso[:,0]==0.)&(collectionColourso[:,1]==0.)&(collectionColourso[:,2]==0.5))[0]
            for j in w:
                collectionColourso[j]=num.array([0.05,0.05,0.05,1.0])

            collection=Poly3DCollection(poly3d,linewidths=1.0,facecolors='k',edgecolors=collectionColoursw)
            collection.set_alpha(1.)
            pylCollection1=ax1.add_collection3d(collection)

    pyl.show()
    sys.exit()

if  saveMovie:
    pyl.ion()

    shapeFig=pyl.figure('Shape',figsize=(10,8))
    shapeFig.subplots_adjust(hspace=0,wspace=0)
    #ax1=shapeFig.add_subplot(121,projection='3d')
    ax1=shapeFig.add_subplot(111,projection='3d')
    ax1.set_aspect('equal')
    ax1.set_xlim(-120,120)
    ax1.set_ylim(-120,120)
    ax1.set_zlim(-120,120)
    ax1.set_axis_bgcolor('0.0')
    ax1.view_init(azim=0, elev=0)
    ax1.set_axis_off()
    pyl.title('H$_2$O Absorption',color='w')

    poly3d=reshape(vertices,vertIndices,
    0,0,0,
    0.,0.,0.)[1]

    w=num.where((collectionColoursw[:,0]==1.)&(collectionColoursw[:,1]==1.))[0]
    for j in w:
        collectionColoursw[j]=num.array([0.05,0.05,0.05,1.0])
    w=num.where((collectionColourso[:,0]==0.)&(collectionColourso[:,1]==0.)&(collectionColourso[:,2]==0.5))[0]
    for j in w:
        collectionColourso[j]=num.array([0.05,0.05,0.05,1.0])

    collection=Poly3DCollection(poly3d,linewidths=1.0,facecolors=collectionColoursw,edgecolors=collectionColoursw)
    collection.set_alpha(1.)
    pylCollection1=ax1.add_collection3d(collection)

    #for animations
    import matplotlib.animation as manimation

    FFMpegWriter = manimation.writers['ffmpeg']#ffmpeg_file
    metadata = dict(title='Phoebe Water-ice', artist='Matplotlib',
    comment='Movie support!')
    #30 FPS
    writer = FFMpegWriter(fps=30, metadata=metadata,bitrate=2400)
    azims=num.arange(0,360*4,1)-45.
    elevs=num.concatenate([num.ones(540)*30.,num.linspace(30,-30,360),num.ones(540)*-30.])
    #15 FPS
    #writer = FFMpegWriter(fps=15, metadata=metadata,bitrate=2400)
    #azims=num.arange(0,360*4,2)-45.
    #elevs=num.concatenate([num.ones(270)*30.,num.linspace(30,-30,180),num.ones(270)*-30.])

    ax1.view_init(azim=0,elev=30.)
    ax1.view_init(azim=10,elev=30.)

    print 'Creating an mpeg'
    stepped=0
    with writer.saving(shapeFig,"/data/PhoebeFigures/Phoebe_test.mp4",300): #dpi specified here
        for i in range(len(azims)):
            print i+1,len(azims),azims[i],elevs[i]
            ax1.view_init(azim=azims[i],elev=elevs[i])
            writer.grab_frame()
            stepped+=1



fig1=pyl.figure(3,figsize=fs)
fig1.subplots_adjust(hspace=0,top=0.99,bottom=0.08)

sp1=fig1.add_subplot(311,axisbg='black')
#if showMap: sp1.imshow(newImage,interpolation='nearest',cmap='gray',extent=extent,zorder=1)
sp1.add_collection(lLCollection1)
mappable1 = pyl.cm.ScalarMappable(cmap = cmap)
mappable1.set_array(num.array([cMinw,cMaxw]))
cbar1 =fig1.colorbar(mappable1,fraction=0.01)#, orientation='horizontal')
cbar1.set_label(cbTitlew)
pyl.ylabel('Latitude (deg)')

if showContours:
    contours=pyl.contour(gridLong,gridLat,gridD-106,numContours,colors='w',zorder=5)
    contours.levels = [nf(val) for val in contours.levels ]
    # Label levels with specially formatted floats
    if pyl.rcParams["text.usetex"]:
        fmt = r'%r'
    else:
        fmt = '%r'
    pyl.clabel(contours, contours.levels, inline=True, fmt=fmt, fontsize=14)


sp2=fig1.add_subplot(312,axisbg='black')#,sharex=sp1,sharey=sp1)
#if showMap: sp2.imshow(newImage,interpolation='nearest',cmap='gray',extent=extent,zorder=1)
sp2.add_collection(lLCollection4)
mappable2 = pyl.cm.ScalarMappable(cmap = cmap)
mappable2.set_array(num.array([cMinwb,cMaxwb]))
cbar2 =fig1.colorbar(mappable2,fraction=0.01)#, orientation='horizontal')
cbar2.set_label(cbTitlewb)
pyl.ylabel('Latitude (deg)')

if showContours:
    contours=pyl.contour(gridLong,gridLat,gridD-106,numContours,colors='w',zorder=5)
    contours.levels = [nf(val) for val in contours.levels ]
    # Label levels with specially formatted floats
    if pyl.rcParams["text.usetex"]:
        fmt = r'%r'
    else:
        fmt = '%r'
    pyl.clabel(contours, contours.levels, inline=True, fmt=fmt, fontsize=14)


sp3=fig1.add_subplot(313)#,sharex=sp1,sharey=sp1)
if showMap: sp3.imshow(newImage,interpolation='nearest',cmap='gray',extent=extent,zorder=1,aspect='auto')
#sp3.add_collection(lLCollection5)
mappable3 = pyl.cm.ScalarMappable(cmap = cmap)
mappable3.set_array(num.array([cMinww,cMaxww]))
cbar3 =fig1.colorbar(mappable3,fraction=0.01)#, orientation='horizontal')
cbar3.set_label(cbTitleww)

#if showContours:
#    contours=pyl.contour(gridLong,gridLat,gridD-106,numContours,colors='w',zorder=5)
#    contours.levels = [nf(val) for val in contours.levels ]
#    # Label levels with specially formatted floats
#    if pyl.rcParams["text.usetex"]:
#        fmt = r'%r'
#    else:
#        fmt = '%r'
#    pyl.clabel(contours, contours.levels, inline=True, fmt=fmt, fontsize=14)


if showJason:
    sp1.set_xlim(-20,80)
    sp1.set_ylim(-20,60)
    sp2.set_xlim(-20,80)
    sp2.set_ylim(-20,60)
    sp3.set_xlim(-20,80)
    sp3.set_ylim(-20,60)
else:
    sp1.set_xlim(-180,180)
    #sp1.set_ylim(-56.9534,56.9534)
    sp1.set_ylim(-90,90)
    sp2.set_xlim(-180,180)
    #sp2.set_ylim(-56.9534,56.9534)
    sp2.set_ylim(-90,90)
    sp3.set_xlim(-180,180)
    #sp3.set_ylim(-56.9534,56.9534)
    sp3.set_ylim(-90,90)
pyl.ylabel('Latitude (deg)')
pyl.xlabel('Longitude (deg)')

if saveMap: pyl.savefig('Phoebe_fullmap.tiff',bbox='tight')
pyl.show()
