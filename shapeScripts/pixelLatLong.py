#! /usr/bin/env python

import commands,sys,jdcal
from astropy.io import fits as pyf

def cutter(s):
    x=s.split()
    try:
        return float(x[2])
    except:
        return x[2]

def bodyCutter(s1,s2):
    x=s1.split()
    y=s2.split()
    return [float(x[2][1:len(x[2])-2]),float(x[3][:len(x[3])-1]),float(y[0][:len(y[0])-1])]


def MJD(x):
    s=x.split('T')
    o=s[0].split('-')
    e=s[1].split(':')

    year=jdcal.gcal2jd(int(float(o[0])),int(float(o[1])),int(float(o[2])))[1]
    year+=float(e[0])/24.+float(e[1])/24./60.+float(e[2])/24./3600.

    return year


def pixLatLong(image):

    with pyf.open(image+'.fits',ignore_missing_end=True) as han:
        data=han[0].data

    (C,A,B)=data.shape

    outObj=[]
    for i in range(A):
        outObj.append([])
        for j in range(B):
            comm='campt from=%s line=%s sample=%s'%(image+'.cub',i+1,j+1)

            output=commands.getoutput(comm).split('\n')
            #dic={'PlanetocentricLatitude':None,
            #     'PositiveEast360Longitude':None,
            #     'SampleResolution':None,
            #     'SubSpacecraftLatitude': None,
            #     'SubSpacecraftLongitude': None,
            #     'SpacecraftAzimuth': None,
            #     'SubSolarLatitude': None,
            #     'SubSolarLongitude': None,
            #     'SubSolarAzimuth': None,
            #     'SubSolarGroundAzimuth': None,
            #     'BodyFixedCoordinate':None,
            #     'SunPosition':None,
            #     'SpacecraftPosition':None}


            dic={'Filename':None,
                'Sample ':None,
                'Line ':None,
                'PixelValue':None,
                'RightAscension':None,
                'Declination':None,
                'PlanetocentricLatitude':None,
                'PlanetographicLatitude':None,
                'PositiveEast360Longitude':None,
                'PositiveEast180Longitude':None,
                'PositiveWest360Longitude':None,
                'PositiveWest180Longitude':None,
                'BodyFixedCoordinate':None,
                'LocalRadius':None,
                'SampleResolution':None,
                'LineResolution':None,
                'SpacecraftPosition':None,
                'SpacecraftAzimuth':None,
                'SlantDistance':None,
                'TargetCenterDistance':None,
                'SubSpacecraftLatitude':None,
                'SubSpacecraftLongitude':None,
                'SpacecraftAltitude':None,
                'OffNadirAngle':None,
                'SubSpacecraftGroundAzimuth':None,
                'SunPosition':None,
                'SubSolarAzimuth':None,
                'SolarDistance':None,
                'SubSolarLatitude':None,
                'SubSolarLongitude':None,
                'SubSolarGroundAzimuth':None,
                'Phase':None,
                'Incidence':None,
                'Emission':None,
                'NorthAzimuth':None,
                'EphemerisTime':None,
                'UTC':None,
                'LocalSolarTime':None,
                'SolarLongitude':None}
                
             
            if '**ERROR**' not in output[0]:
                print i,j
                for k in range(len(output)):
                    print output[k]
                    for l in dic:
                        if l in output[k]:
                            if 'BodyFixedCoordinate' in l or 'SunPosition' in l or 'SpacecraftPosition' in l:
                                dic[l]=bodyCutter(output[k],output[k+1])
                            else:   
                                dic[l]=cutter(output[k])
                print '\n\n\n'
            #print dic
            outObj[i].append(dic.copy())
    return outObj

if __name__=="__main__":
    import pickle as pick,pylab as pyl,numpy as num,os
    from lineFitUtils import poly
    fn='cv1465671593_1'
    image=fn+'_ir'
    dir='2004163T121836_2004163T192848'
    #os.chdir(dir)
    out=pixLatLong(image)
    with open(image+'.campt.pickle','w+') as han:
        pick.dump(out,han)
        
    image=fn+'_vis'
    out=pixLatLong(image)
    with open(image+'.campt.pickle','w+') as han:
        pick.dump(out,han)

    #os.chdir('..')
    sys.exit()

    
    with open(image+'.campt.pickle') as han:
        out=pick.load(han)
        

    Cass=[]
    Phoebe=[]
    for i in range(len(out)):
        for j in range(len(out[i])):
            if out[i][j]['SpacecraftPosition']<>None:
                [X,Y,Z]=out[i][j]['SpacecraftPosition']
                Cass.append([X,Y,Z,MJD(out[i][j]['UTC']),out[i][j]['TargetCenterDistance'],out[i][j]['SubSolarLongitude'],out[i][j]['SubSolarLatitude'],out[i][j]['SubSolarAzimuth'],i,j])
                [X,Y,Z]=out[i][j]['BodyFixedCoordinate']
                Phoebe.append([X,Y,Z])
                print out[i][j]['SubSolarLongitude'],out[i][j]['SubSolarLatitude'],out[i][j]['SubSolarAzimuth']
    Cass=num.array(Cass)
    args=num.argsort(Cass[:,3])
    Cass=Cass[args]
    Phoebe=num.array(Phoebe)[args]
    #print Cass[0]
    #print Cass[len(Cass)-1]
    #sys.exit()
    
    order=2
    objx=poly(Cass[:,3],Cass[:,0],renormX=True,order=order)
    objy=poly(Cass[:,3],Cass[:,1],renormX=True,order=order)
    objz=poly(Cass[:,3],Cass[:,2],renormX=True,order=order)

    
    fig1=pyl.figure(1)

    pyl.plot(Cass[:,9],Cass[:,7])
    pyl.show()

