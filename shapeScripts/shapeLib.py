import numpy as num,sys
from numpy import linalg
from timeit import default_timer as timer
from numbapro import vectorize,guvectorize
from math import cos,acos

d2r=num.pi/180.
r2d=180./num.pi
twopi=2*num.pi

def mag(inVec):
    return num.sum(inVec**2)**0.5

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

def rotx(ang):
        
    w=num.array([[1.0,0.0,0.0],
                [0.0,num.cos(ang),-num.sin(ang)],
                [0.0,num.sin(ang),num.cos(ang)]])
    return w
def roty(ang):
    w=num.array([[num.cos(ang),0.0,-num.sin(ang)],
                [0.0,1.0,0.0],
                [num.sin(ang),0.0,num.cos(ang)]])
    return w

def rotz(ang):
    w=num.array([[num.cos(ang),-num.sin(ang),0.0],
                [num.sin(ang),num.cos(ang),0.0],
                [0.0,0.0,1.0]])
    return w

@vectorize(['float64(float64, float64, float64)'],target='cpu')
def magVectorize(a,b,c):
    return (a**2+b**2+c**2)**0.5
    
def arrayMagRepVectorize(inVec):
    (AA,BB,CC)=inVec.shape
    out=num.zeros((AA,BB)).astype(inVec.dtype)
    out[:,0]=magVectorize(inVec[:,0,0],inVec[:,0,1],inVec[:,0,2])
    out[:,1]=magVectorize(inVec[:,1,0],inVec[:,1,1],inVec[:,1,2])
    out[:,2]=magVectorize(inVec[:,2,0],inVec[:,2,1],inVec[:,2,2])
    return num.repeat(out,CC,axis=1).reshape(AA,BB,CC)

def arrayMagRep(inVec): 
    (AA,BB,CC)=inVec.shape
    x=inVec**2
    mags=num.sum(x,axis=2)**0.5
    repMags=num.repeat(mags,CC,axis=1).reshape(AA,BB,CC)
    return repMags

def arrayMag2Vectorize(inVec):
    (AA,BB)=inVec.shape
    out=num.zeros(AA).astype(inVec.dtype)
    out[:]=magVectorize(inVec[:,0],inVec[:,1],inVec[:,2])
    return out

def arrayMag2(inVec):
    (AA,BB)=inVec.shape
    x=inVec**2
    mags=num.sum(x,axis=1)**0.5
    return mags

@vectorize(['float64(float64, float64, float64)'], target='cpu')
def getM(a,b,c):
    return (a+b+c)/3.
def getMids(p):
    (AA,BB,CC)=p.shape
    out=num.zeros((AA,BB)).astype(p.dtype)
    out[:,0]=getM(p[:,0,0],p[:,1,0],p[:,2,0])
    out[:,1]=getM(p[:,0,1],p[:,1,1],p[:,2,1])
    out[:,2]=getM(p[:,0,2],p[:,1,2],p[:,2,2])
    return out

"""
@vectorize(['float64(float64, float64,float64, float64)'],target='cpu')
def cr0(a,b,c,d):
    return a*d-b*c
@vectorize(['float64(float64, float64,float64, float64)'],target='cpu')
def cr1(a,b,c,d):
    return -a*d+b*c
def crossVectorize(u,v):
    out=num.zeros(u.shape).astype(u.dtype)

    out[:,0]=cr0(u[:,1],u[:,2],v[:,1],v[:,2])
    out[:,1]=cr1(u[:,0],u[:,2],v[:,0],v[:,2])
    out[:,2]=cr0(u[:,0],u[:,1],v[:,0],v[:,1])
    return out
"""
def midsNormals(poly3d):
    mids=num.array(getMids(poly3d))
    #mids=num.mean(poly3d,axis=1) #midpoints
    (A,B)=mids.shape
    repMids=num.repeat(mids,3,axis=0).reshape(A,B,3) 
    abc=poly3d-repMids  #non-normalized abc vectors
    #abc/=arrayMagRep(abc)  #abc now normzalized
    abc/=num.array(arrayMagRepVectorize(abc))  #abc now normzalized
    a=abc[:,0]
    b=abc[:,1]
    c=abc[:,2]

    ###corss products take longer with vectorize
    axb=num.cross(a,b)
    bxc=num.cross(b,c)
    cxa=num.cross(c,a)
    cross=num.concatenate([num.array([axb]).transpose(),num.array([bxc]).transpose(),num.array([cxa]).transpose()]).transpose().reshape(A,3,3) #unnormalized cross product vectors
    #cross/=arrayMagRep(cross) #normalized cross products
    cross/=num.array(arrayMagRepVectorize(cross)) #normalized cross products
    
    normals=num.mean(cross,axis=1) #normals

    return (mids,normals)


            
@vectorize(['float64(float64, float64, float64, float64, float64, float64)'], target='cpu')
def dot1(w0,w1,w2,v0,v1,v2):
    return w0*v0+w1*v1+w2*v2

def rotDot(W,vertices):
    out=num.zeros(vertices.shape).astype(vertices.dtype)
    out[:,0]=dot1(W[0][0],W[0][1],W[0][2],vertices[:,0],vertices[:,1],vertices[:,2])
    out[:,1]=dot1(W[1][0],W[1][1],W[1][2],vertices[:,0],vertices[:,1],vertices[:,2])
    out[:,2]=dot1(W[2][0],W[2][1],W[2][2],vertices[:,0],vertices[:,1],vertices[:,2])
    return out

def reshape(vertices,vertIndices,
            phi_o=200.35216011092,theta_o=-21.76086293882,az_o=264.28748413648,
            phi_s=291.5572425051,theta_s=-12.85085685023,az_s=90.17):
            #offsets=num.array([0.,0.])):
    
    """
            Pass the actual camInfo log values (in degrees and with the sign in the log file!)

            
    defaults are from ISS image cN1465652275_1.caminfo
    ####it seems we need negatives infront of all the angles!
    longitude phi_o=d2r*-(200.35216011092)
    latitude theta_o=d2r*-(-21.76086293882)
    azimuth az_o=d2r*-( 264.28748413648)
    """


    
    Aphi_s=phi_s*-d2r
    Aphi_o=phi_o*-d2r
    Atheta_s=theta_s*-d2r
    Atheta_o=theta_o*-d2r
    #azimuthes rotate along vector out of the screen
    ###this current orientation matches the ds9 image when x is swapped
    Aaz_o=az_o*-d2r
    Aaz_s=az_s*-d2r

    n_obs_r=num.array([num.cos(Aphi_o)*num.cos(Atheta_o),
                      num.sin(Aphi_o)*num.cos(Atheta_o),
                      num.sin(Atheta_o)])
    n_sun=num.array([num.cos(Aphi_s)*num.cos(Atheta_s),
                    num.sin(Aphi_s)*num.cos(Atheta_s),
                    num.sin(Atheta_s)])*1.4e9

    n_obs=num.array([1.,0.,0.])*1.4e9


    #apply a triplet of rotations to the body to get the geometry correct
    #longitude pivots around the z axis and has zero at +i (xhat)
    #latitude pivots around the y axis and has zero in the xy plane
    #azimuth is about the z axis again

    w1=rotz(Aphi_o)
    w2=roty(Atheta_o)
    w3=rotx(Aaz_o)
    #w1=rot(Aphi_o,'z')
    #w2=rot(Atheta_o,'y')
    #w3=rot(Aaz_o,'x')
    W=num.dot(w3,num.dot(w2,w1))


    w3=rotx((Aaz_s-Aaz_o))
    rot_n_sun=num.dot(w3, num.array([num.cos((Aphi_s-Aphi_o))*num.cos((Atheta_s-Atheta_o)),
                                     num.sin((Aphi_s-Aphi_o))*num.cos((Atheta_s-Atheta_o)),
                                     num.sin((Atheta_s-Atheta_o))])*1.49e9)
    print rot_n_sun,'*'

    rot_vertices=num.array(rotDot(W,vertices))

    #vertices/rot_vertices is just a list of vertices while vertIndices is a list of triplets
    #of the three vertices that define a facet.

    poly3d=num.copy(rot_vertices[vertIndices])

    (mids,normals)=midsNormals(poly3d)
    #print mids,normals.shape

    #w=num.where(arrayMag2(mids-normals)>arrayMag2(mids+normals))
    w=num.where(arrayMag2Vectorize(mids-normals)>arrayMag2Vectorize(mids+normals))
    normals[w]*=-1
    #### end replacement

    
    #get the observability
    #illuminated, visible, and both facets
    mag_n_sun=num.sum(rot_n_sun**2)**0.5
    mag_n_obs=num.sum(n_obs**2)**0.5

    #print normals.shape,rot_n_sun.shape
    #angs_s=num.arccos(num.dot(normals,rot_n_sun)/mag_n_sun)
    angs_s=num.array(arccosDot(normals[:,0],normals[:,1],normals[:,2],rot_n_sun[0],rot_n_sun[1],rot_n_sun[2],mag_n_sun))
    illuminated=num.where((angs_s<=num.pi/2)&(angs_s>=-num.pi/2))[0]

    #angs_o=num.arccos(num.dot(normals,n_obs)/mag_n_obs)
    angs_o=num.array(arccosDot(normals[:,0],normals[:,1],normals[:,2],n_obs[0],n_obs[1],n_obs[2],mag_n_obs))
    #visible=num.where((angs_o>=-num.pi/2)&(angs_o<=num.pi/2))[0]

    ill_and_obs=num.where( ((angs_s>=-num.pi/2)&(angs_s<=num.pi/2))
                        & ((angs_o>=-num.pi/2)&(angs_o<=num.pi/2)))[0]




    #setup the colours array
    #0.0 is black,1.0 is white
    ###since the fitting code doesn't want colour but just illuminated or not, we could just do colours=1 where ill_and_obs==True, may be a little faster
    #colours=num.repeat(num.cos(angs_o)*num.cos(angs_s),3).reshape(len(mids),3)


    
    return (rot_vertices,num.copy(poly3d),num.copy(normals),num.copy(mids),num.copy(ill_and_obs),n_obs,rot_n_sun,angs_s,angs_o)
    #(colours,area)=getColourArea(angs_o,angs_s,mids,poly3d)
    #return (rot_vertices,num.copy(poly3d),num.copy(colours),num.copy(area),num.copy(normals),num.copy(mids),num.copy(ill_and_obs),n_obs,rot_n_sun,angs_s,angs_o)

def getColourArea(angs_o,angs_s,mids,poly3d):
    colours=num.repeat(getc(angs_o,angs_s),3).reshape(len(mids),3)
    num.clip(colours,0.0,num.max(colours),colours)
    not_ill_and_obs=num.where(((angs_s>=num.pi/2)|(angs_s<=-num.pi/2))
                            &((angs_o>=num.pi/2)|(angs_o<=-num.pi/2)))[0]   ####and and or produce the same illuminated area
    colours[not_ill_and_obs]=0.0
    #below two lines are a visualization hack
    #ill_and_obs=num.where( ((angs_s>=-num.pi/2)&(angs_s<=num.pi/2))
    #                    & ((angs_o>=-num.pi/2)&(angs_o<=num.pi/2)))[0]
    #colours[ill_and_obs]=1.0
    cmax=num.max(colours)
    colours/=cmax
    ####


    #l1=((poly3d[:,0,0]-poly3d[:,1,0])**2+(poly3d[:,0,1]-poly3d[:,1,1])**2)**0.5
    #l2=((poly3d[:,0,0]-poly3d[:,2,0])**2+(poly3d[:,0,1]-poly3d[:,2,1])**2)**0.5
    #l3=((poly3d[:,1,0]-poly3d[:,2,0])**2+(poly3d[:,1,1]-poly3d[:,2,1])**2)**0.5
    #s=(l1+l2+l3)/2.
    #area=(s*(s-l1)*(s-l2)*(s-l3))**0.5
    l1=getL(poly3d[:,0,0],poly3d[:,1,0],poly3d[:,0,1],poly3d[:,1,1])
    l2=getL(poly3d[:,0,0],poly3d[:,2,0],poly3d[:,0,1],poly3d[:,2,1])
    l3=getL(poly3d[:,1,0],poly3d[:,2,0],poly3d[:,1,1],poly3d[:,2,1])
    s=lSum(l1,l2,l3)
    area=getArea(s,l1,l2,l3)
    area[num.where(num.isnan(area))]=0.
    return (colours,area)

@vectorize (['float64(float64, float64, float64, float64, float64, float64, float64)'],target='cpu')
def arccosDot(a,b,c,d,e,f,mag):
    return acos((a*d+b*e+c*f)/mag)

@vectorize (['float64(float64, float64)'],target='cpu')
def getc(a_o,a_s):
    return cos(a_o)*cos(a_s)

@vectorize (['float64(float64, float64, float64, float64)'], target='cpu')
def getL(a,b,c,d):
    return ((a-b)*(a-b)+(c-d)*(c-d))**0.5

@vectorize (['float64(float64, float64, float64)'], target='cpu')
def lSum(l1,l2,l3):
    return (l1+l2+l3)/2.

@vectorize (['float64(float64, float64, float64, float64)'], target='cpu')
def getArea(s,l1,l2,l3):
    return (s*(s-l1)*(s-l2)*(s-l3))**0.5

def shapeGen_ISS(vertices,vertIndices,
                long_o,lat_o,az_o,
                long_s,lat_s,az_s,
                offsets,sampleResolution,
                imData):

    x=reshape(vertices,vertIndices,
                phi_o=long_o,theta_o=lat_o,az_o=az_o,
                phi_s=long_s,theta_s=lat_s,az_s=az_s)
                #offsets=offsets)

   
    (rot_vertices,poly3d,n,m,ill_and_obs,n_obs,rot_n_sun,angs_s,angs_o)=x
    (colours,area)=getColourArea(angs_o,angs_s,m,poly3d)

    p3d=poly3d[ill_and_obs]
    c3d=colours[ill_and_obs]
    c2d=num.copy(c3d)
    m2d=m[ill_and_obs]
    a2d=area[ill_and_obs]

    #(rot_vertices,poly3d,c3d,m2d,c2d,a2d,n,m,n_obs,rot_n_sun,angs_s,angs_o)=x
    #return (rot_vertices,num.copy(poly3d[ill_and_obs]),num.copy(colours[ill_and_obs]),num.copy(m[ill_and_obs]),colours2d,area,n,m,n_obs,rot_n_sun,angs_s,angs_o)

    #generate the model image
    c2d=c2d[:,0]
    c2d*=a2d
    c2d/=num.max(c2d)

    (A,B)=imData.shape
    y=num.arange(-A*sampleResolution/2,A*sampleResolution/2,sampleResolution)
    z=num.arange(-B*sampleResolution/2,B*sampleResolution/2,sampleResolution)
    image=num.zeros([len(y),len(z)]).astype('float')
    
    K=((m2d[:,1]-y[0])/(y[1]-y[0])).astype('int')
    L=((m2d[:,2]-z[0])/(z[1]-z[0])).astype('int')

    K+=int(offsets[0])
    L+=int(offsets[1])
    
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





def shapeGen_VIMS(vertices,vertIndices,
                  long_o_not,lat_o_not,az_o_not,
                  long_s,lat_s,az_s,
                  distancenot,
                  offsets, offsetVel,offsetVelAngle,
                  pixScale,
                  inot,jnot,
                  deltaX,deltaY,deltaZ,
                  lons,lats,
                  pixelTimes,
                  imData,verbose=False,vis=False,mask=None,
                  az_adjust=0.0):

    w1=rot(long_o_not*-d2r,'z')
    w2=rot(lat_o_not*-d2r,'y')
    w3=rot((az_o_not+az_adjust)*-d2r,'x')
    W_not=num.dot(w3,num.dot(w2,w1))


    #need to rotate the spaceCraftVector according to the free params,
    #subspacecraft long,lat,az because the one in the image header
    #are incorrect once these values are chosen
    W_not_reverse=linalg.inv(W_not)
    spaceCraftVectornot=num.dot(W_not_reverse,distancenot*num.array([1.,0.,0.]))
    
    x=reshape(vertices,vertIndices,
              long_o_not,lat_o_not,az_o_not+az_adjust,
              long_s,lat_s,az_s+az_adjust)

    (rot_vertices_not,poly3d,n_not,m_not,ill_and_obs,n_obs,rot_n_sun,angs_s,angs_o)=x
    #(colours,area)=getColourArea(angs_o,angs_s,m_not,poly3d)

    normMids_not=m_not/num.repeat((m_not[:,0]**2+m_not[:,1]**2+m_not[:,2]**2)**0.5,3).reshape(m_not.shape[0],3)
    if verbose: print rot_vertices_not[0],long_o_not,lat_o_not
    
    velAng=num.array([num.cos(offsetVelAngle*d2r),num.sin(offsetVelAngle*d2r)])

    (A,B)=imData.shape
    image=imData*0.0


    sampleResolutionnot=pixScale*distancenot/1000.


    #rotated to have the planet vertex at 0.,0.
    #now need to get the angle between the next i,j combination and the spaceCraftVectornot
    #the next i,j space craft vector is spaceCraftVectornot+deltaXYZ[i,j]
    currLon,currLat,currRes=long_o_not,lat_o_not,sampleResolutionnot
    vertsInPixel=[]
    for i in range(A): #line
        vertsInPixel.append([])
        for j in range(B): #sample

            #rotate the current spaceCraftVector by the not angles
            #this places the planet vertex at 0.0
            spaceCraftVector=num.dot(W_not, spaceCraftVectornot+num.array([deltaX[i,j],deltaY[i,j],deltaZ[i,j]]))
            spaceCraftDistance=mag(spaceCraftVector)
            spaceCraftVector/=spaceCraftDistance
            
            #d=num.dot(normMids_not,spaceCraftVector)
            d=dot1(normMids_not[:,0],normMids_not[:,1],normMids_not[:,2],spaceCraftVector[0],spaceCraftVector[1],spaceCraftVector[2])
            planetVertex=num.argmax(d)

            sampleResolution=pixScale*spaceCraftDistance/1000.
            
        
        
            #print i,j,sampleResolution,spaceCraftDistance,
            #the resolution change is rarely the issue
            if lons[planetVertex]<>currLon or lats[planetVertex]<>currLat or num.abs(sampleResolution-currRes)/currRes>0.01:
                currLon,currLat,currRes=lons[planetVertex],lats[planetVertex],sampleResolution
                x=reshape(vertices,vertIndices,
                          currLon,currLat,az_o_not+az_adjust,
                          long_s,lat_s,az_s+az_adjust)
                (rot_vertices,poly3d,n,m,ill_and_obs,n_obs,rot_n_sun,angs_s,angs_o)=x
                 
                if verbose: print rot_vertices[0],currLon,currLat,long_o_not,lat_o_not
                #now we have the shape model rotated such that the desired subspacecraft lat/long are aligned with (0.,0.,1.)
                #next we spatially offset the model in the x/y plane by offsets (x,y) where we always assume that pixel 0,0 is offset 0,0
                #then we image the pixel i,j

            (offy,offz)=offsets+offsetVel*(pixelTimes[i,j]-pixelTimes[inot,jnot])*24.*3600.*velAng


            #pixeli,j are the distance from the reference point on the body in the image plane
            if vis:
                sampResHor=sampleResolution
                sampResVert=sampleResolution
            else:
                sampResHor=sampleResolution
                sampResVert=sampleResolution/2.
            
            pixel_off_i=(i-0)*sampResHor
            pixel_off_j=(j-0)*sampResVert
            #print pixel_off_i,pixel_off_j

        
            m2d=num.copy(m)[ill_and_obs]


            #offset by the difference between position offset and the pixel offset
            w_in_pixel=num.where((num.abs(m2d[:,1]+(offy-pixel_off_i))<(sampResHor/2.)) & (num.abs(m2d[:,2]+(offz-pixel_off_j))<(sampResVert/2.)))
            #print len(w_in_pixel[0])
            if len(w_in_pixel[0])>0:
                image[i,j]=1.0
                vertsInPixel[i].append(ill_and_obs[w_in_pixel]*1.0)
            else:
                vertsInPixel[i].append(num.array([]))
                

            if i==inot and j==jnot:
                poly3d_out=num.copy(poly3d)
                (colours_out,area)=getColourArea(angs_o,angs_s,m,poly3d_out)#=num.copy(colours)
                rot_vert=num.copy(rot_vertices)
    image=image[::-1,:]

    #edgeChi(imData,image)

    #from the cub image (imData), 256 is illuminated, 0 is dark
    #from the model image (image), 1 is illuminated, 0 is dark
    if vis:
        chi=-num.sum(((imData/256.)-(image*mask))**2)#/(A*B-7.)
    else:
        chi=-num.sum(((imData/256.)-(image))**2)#/(A*B-7.)
    #print chi
    #chi=edgeChi(imData,image)
    if verbose: print chi
    return (image,poly3d_out,colours_out,rot_vert,vertsInPixel,chi)
    

def edgeChi(im1,im2):
    i1=im1/num.max(im1)
    i2=im2/num.max(im2)

    im1pix=[]
    for ii in range(len(i1)):
        w=num.where(i1[ii]>0) #all illuminated pixels in this line
        #print ii,w
        if len(w[0])>0:
            W=num.sort(w[0])
            if len(W)<4:
                for jj in range(len(W)):
                    im1pix.append([ii,W[jj]])
            else:
                im1pix.append([ii,W[0]])
                im1pix.append([ii,W[1]])
                im1pix.append([ii,W[len(W)-2]])
                im1pix.append([ii,W[len(W)-1]])


    im2pix=[]
    for ii in range(len(i1)):
        w=num.where(i2[ii]>0) #all illuminated pixels in this line
        #print ii,w
        if len(w[0])>0:
            W=num.sort(w[0])
            if len(W)<4:
                for jj in range(len(W)):
                    im2pix.append([ii,W[jj]])
            else:
                im2pix.append([ii,W[0]])
                im2pix.append([ii,W[1]])
                im2pix.append([ii,W[len(W)-2]])
                im2pix.append([ii,W[len(W)-1]])

    I1=im1*0.0
    for ii in range(len(im1pix)):
        I1[im1pix[ii][0],im1pix[ii][1]]=1.0
    I2=im2*0.0
    for ii in range(len(im2pix)):
        I2[im2pix[ii][0],im2pix[ii][1]]=1.0
    #print I1-I2
    chi=-num.sum(((I1)-(I2))**2)
    return chi
    #print chi,'**'
    #print num.sum((i1-i2)**2)
    #sys.exit()
