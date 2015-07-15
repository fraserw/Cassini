#! /usr/bin/env python

import glob
import commands,os
import sys
import numpy as num, scipy as sci,pickle as pick
from astropy.io import fits as pyf


def cmp0(a,b):
    return cmp(a[0],b[0])

def findPhoebeCubes(dir):
    cubes=[]    
    files=glob.glob(dir+'/*lbl')
    for ii in files:
        with open(ii) as han:
            data=han.readlines()
            for jj in range(len(data)):
                if 'TARGET_NAME' in data[jj]:
                    target=data[jj].split('=')[1].split()[0][1:]
                    target=target[:len(target)-1]
                    if target=='PHOEBE':
                        cubes.append(ii.split('.')[0])
                    break
    if cubes==[]:
        return None
    else:
        return cubes

def preprocCubes(cubeName,procDir='/data/VIMS/covims_0004/procdata'):
    if not os.path.exists(procDir):
        os.system('mkdir %s'%(procDir))

    cwd=os.getcwd()

    cube=cubeName.split('/')[1].split('.')[0]
    
    opDir=procDir+'/'+cubeName.split('/')[0]
    if not os.path.exists(opDir):
        os.system('mkdir %s'%(opDir))

    cpComm='cp %s* %s'%(cubeName,opDir)
    print cpComm
    os.system(cpComm)

    os.chdir(opDir)

    convComm='vims2isis from=%s.qub vis=%s_vis.cub ir=%s_ir.cub'%(cube,cube,cube)
    print
    print ' >%s'%(convComm)
    convCommOut=commands.getoutput(convComm)
    print convCommOut
    
    specCommVis='spiceinit from=%s_vis.cub'%(cube)
    specCommIR ='spiceinit from=%s_ir.cub'%(cube)
    specCommVisOut=commands.getoutput(specCommVis)
    specCommIROut =commands.getoutput(specCommIR)

    print
    print ' >%s'%(specCommVis)
    print specCommVisOut
    print
    print ' >%s'%(specCommIR)
    print specCommIROut

    vimscalVisComm='vimscal from=%s_vis.cub to=c%s_vis.cub'%(cube,cube)
    vimscalIRComm ='vimscal from=%s_ir.cub to=c%s_ir.cub'%(cube,cube)
    vimscalVisOut=commands.getoutput(vimscalVisComm)
    vimscalIROut =commands.getoutput(vimscalIRComm)

    print
    print ' >%s'%(vimscalVisComm)
    print vimscalVisOut
    print
    print ' >%s'%(vimscalIRComm)
    print vimscalIROut

    try:
        os.remove('%s.qub'%cube)
    except:
        pass
    
    os.chdir(cwd)

    log={}
    log[0]=convCommOut
    log[1]=specCommVisOut
    log[2]=specCommIROut
    log[3]=vimscalVisOut
    log[4]=vimscalIROut

    return log


def genFitsMeans(cubeName,procDir='/data/VIMS/covims_0004/procdata'):
    log={}
    cwd=os.getcwd()

    cube=cubeName.split('/')[1].split('.')[0]
    opDir=procDir+'/'+cubeName.split('/')[0]

    os.chdir(opDir)

    try:
        with pyf.open('c'+cube+'_ir.fits') as irHan:
            ir=irHan[0].data
            mean=num.mean(num.concatenate([ir[:34,:,:],ir[82:105]]),axis=0)
            print mean.shape
            HDU=pyf.PrimaryHDU(mean)
            list=pyf.HDUList([HDU])
            if os.path.isfile('c'+cube+'_ir_mean.fits'):
                os.remove('c'+cube+'_ir_mean.fits')
            list.writeto('c'+cube+'_ir_mean.fits')
    except: pass
    try:
        with pyf.open('c'+cube+'_vis.fits') as visHan:
            vis=visHan[0].data
            mean=num.mean(vis[:100,:,:],axis=0)
            HDU=pyf.PrimaryHDU(mean)
            list=pyf.HDUList([HDU])
            if os.path.isfile('c'+cube+'_vis_mean.fits'):
                os.remove('c'+cube+'_vis_mean.fits')
            list.writeto('c'+cube+'_vis_mean.fits')
    except: pass
                         


    os.chdir(cwd)

    
def convert2Fits(cubeName,procDir='/data/VIMS/covims_0004/procdata'):
    log={}
    cwd=os.getcwd()

    cube=cubeName.split('/')[1].split('.')[0]
    opDir=procDir+'/'+cubeName.split('/')[0]

    os.chdir(opDir)

    convCommIR='isis2fits from=c%s_ir.cub to=c%s_ir.fits'%(cube,cube)
    print
    print ' >%s '%(convCommIR)
    out=commands.getoutput(convCommIR)
    print out
    print
    log[0]=out

    convCommVis='isis2fits from=c%s_vis.cub to=c%s_vis.fits'%(cube,cube)
    print
    print ' >%s '%(convCommVis)
    out=commands.getoutput(convCommVis)
    print out
    print
    log[1]=out
    
    os.chdir(cwd)
    return log

def getCamStats(cubeName,procDir='/data/VIMS/covims_0004/procdata'):
    log={}
    cwd=os.getcwd()

    cube=cubeName.split('/')[1].split('.')[0]
    opDir=procDir+'/'+cubeName.split('/')[0]

    os.chdir(opDir)

    log[0]=cubeName

    camStatCommVis='camstats from=%s_vis.cub'%(cube)
    camStatCommIR='camstats from=%s_ir.cub'%(cube)

    outVis=commands.getoutput(camStatCommVis)
    outIR=commands.getoutput(camStatCommIR)

    log[1]=outVis
    log[2]=outIR

    os.chdir(cwd)

    return log

def getCamInfo(cubeName,procDir='/data/VIMS/covims_0004/procdata'):
    log={}
    cwd=os.getcwd()

    cube=cubeName.split('/')[1].split('.')[0]
    opDir=procDir+'/'+cubeName.split('/')[0]

    os.chdir(opDir)

    log[0]=cubeName

    
    camStatCommVis='caminfo from=%s_vis.cub camstats=yes to=%s_vis.caminfo'%(cube,cube)
    camStatCommIR='caminfo from=%s_ir.cub camstats=yes to=%s_ir.caminfo'%(cube,cube)

    #print os.system('ls *caminfo')

    
    outVis=[]
    outIR=[]

    outVis=commands.getoutput(camStatCommVis)
    outIR=commands.getoutput(camStatCommIR)

    with open('%s_vis.caminfo'%(cube)) as o:
        log[1]=o.readlines()
    with open('%s_ir.caminfo'%(cube)) as o:
        log[2]=o.readlines()


    os.chdir(cwd)

    return log


def parseLog(x):
    mainKeys=[]
    mainDict={}
    for ii in range(len(x)):
        if x[ii][:11]=='  Object = ':
            s=x[ii].split()
            mainKeys.append(s[2])
            dict={}
            for jj in range(ii,len(x)):
                s=x[jj].split()
                key=s[0]
                try: val=float(s[2])
                except: val=None

                dict[key]=val
                if 'End_Object' in x[jj]: break
            mainDict[mainKeys[len(mainKeys)-1]]=dict.copy()
    return mainDict
    

def parseCamStatsLog():
    camStatsIR=[]
    camStatsVis=[]
    with open('/data/VIMS/covims_0004/procdata/camStats.log') as logHan:
        log=logHan.readlines()
        for i in range(len(log)):
            if '/' in log[i] and '_' in log[i]:
                s=log[i].split('/')
                dir=s[0]
                fn='c'+s[1][:len(s[1])-1]

                imageVis={'dir':dir, 'fileName':fn+'_vis', 'LatitudeAverage':None, 'LongitudeAverage':None, 'ResolutionAverage':None, 'PhaseAverage':None, 'EmissionAverage':None, 'IncidenceAverage':None}
                for j in range(1,91):
                    if '  ' in log[i+j]:
                        s=log[i+j].split()
                        for k in imageVis:
                            if s[0]==k:
                                if s[2]<>'NULL':
                                    imageVis[k]=float(s[2])


                imageIR={'dir':dir, 'fileName':fn+'_ir', 'LatitudeAverage':None, 'LongitudeAverage':None, 'ResolutionAverage':None, 'PhaseAverage':None, 'EmissionAverage':None, 'IncidenceAverage':None}
                for j in range(92,92+90):
                    if '  ' in log[i+j]:
                        s=log[i+j].split()
                        for k in imageIR:
                            if s[0]==k:
                                if s[2]<>'NULL':
                                    imageIR[k]=float(s[2])
                camStatsIR.append(imageIR)
                camStatsVis.append(imageVis)
    return (camStatsVis,camStatsIR)



def parseCamInfoLog():
    """
    This is the more detailed log containing the necessary camera and Sun geometries.
    """
    
    camInfoIR=[]
    camInfoVis=[]
    with open('/data/VIMS/covims_0004/procdata/camInfo.log') as logHan:
        log=logHan.readlines()

    for i in range(len(log)):
        if '/' in log[i] and '_' in log[i]:
            s=log[i].split('/')
            dir=s[0]
            fn='c'+s[1][:len(s[1])-1]

            imageVis={'dir':dir, 'fileName':fn+'_vis', 'SubSolarAzimuth':None, 'SubSolarLatitude':None, 'SubSolarLongitude':None, 'SampleResolution':None, 'SubSpacecraftAzimuth':None, 'SubSpacecraftLongitude':None, 'SubSpacecraftLatitude':None}
            j=1
            while log[i+j]<>'End\n':
                if '  ' in log[i+j]:
                    s=log[i+j].split()
                    for k in imageVis:
                        if s[0]==k:
                            if s[2]<>'NULL':
                                imageVis[k]=float(s[2])
                j+=1

            imageIR={'dir':dir, 'fileName':fn+'_ir', 'SubSolarAzimuth':None, 'SubSolarLatitude':None, 'SubSolarLongitude':None, 'SampleResolution':None, 'SubSpacecraftAzimuth':None, 'SubSpacecraftLongitude':None, 'SubSpacecraftLatitude':None}
            while i+j<len(log)-1 and log[i+j]<>'End\n':
                print len(log),i+j
                if '  ' in log[i+j]:
                    s=log[i+j].split()
                    for k in imageIR:
                        if s[0]==k:
                            if s[2]<>'NULL':
                                imageIR[k]=float(s[2])
                j+=1
                """
                shape=None
                try:
                    imHan=pyf.open(dir+'/'+fn+'_vis_mean.fits')
                    shape=imHan[0].data.shape
                    imHan.close
                except: pass
                imageVis['shape']=shape
                shape=None
                try:
                    imHan=pyf.open(dir+'/'+fn+'_ir_mean.fits')
                    shape=imHan[0].data.shape
                    imHan.close()
                except: pass
                imageIR['shape']=shape
                """
                    
            camInfoIR.append(imageIR)
            camInfoVis.append(imageVis)
            if imageIR['SubSolarAzimuth']<>None:
                print imageIR
                sys.exit()
    return (camInfoVis,camInfoIR)


          
if __name__=="__main__":
    #run this from the covims data directory!
    #/Volumes/data/VIMS/covims_0004/data  

    if os.getcwd()<>'/data/VIMS/covims_0004/data':
        print 'This needs to be run in /data/VIMS/covims_0004/data ONLY!'
        sys.exit()
        
        
    #get Phoebe cubes
    dirs=glob.glob('*')
    cubes=[]
    for i in range(len(dirs)):
        out=findPhoebeCubes(dirs[i])
        if out<>None:
            for jj in range(len(out)):
                cubes.append(out[jj])

    """
    #pre processing all Phoebe data
    cubeProcLog=open('/data/VIMS/covims_0004/procdata/cubeProc.log','w+')    

    for i in range(len(cubes)):
        log=preprocCubes(cubes[i])
        print >>cubeProcLog,cubes[i]
        for j in [0,1,2,3,4]:
            print >>cubeProcLog,log[j]
            print >>cubeProcLog
        print >>cubeProcLog
    cubeProcLog.close()

                
    #convert to fits
    cubeProcLog=open('/data/VIMS/covims_0004/procdata/cubeProc.log','a+')
    print >> cubeProcLog
    for i in range(len(cubes)):
        log=convert2Fits(cubes[i])
        for j in [0,1]:
            print >>cubeProcLog,log[j]
            print >>cubeProcLog
    cubeProcLog.close()

    """
    #gen fits means
    for i in range(len(cubes)):
        genFitsMeans(cubes[i])
    sys.exit()
    """

        

    #####NO LONGER USED
    #create the camstats log
    camStatLog=open('/data/VIMS/covims_0004/procdata/camStats.log','w+')
    for i in range(len(cubes)):
        log=getCamStats(cubes[i])
        print >> camStatLog,log[0]
        print >> camStatLog,log[1]
        print >> camStatLog
        print >> camStatLog,log[2]
        print >> camStatLog,'\n\n'
    camStatLog.close()

    #parse the camera statistics output
    (camStatsVis,camStatsIR)=parseCamStatsLog()
    camStatPickle=open('/data/VIMS/covims_0004/procdata/camStats.pick','w+')
    pick.dump([camStatsVis,camStatsIR],camStatPickle)
    camStatPickle.close()
    ######


    
    
    #create the camInfo log
    camInfoLog=open('/data/VIMS/covims_0004/procdata/camInfo.log','w+')
    visObj=[]
    IRObj=[]
    for i in range(len(cubes)):
        s=cubes[i].split('/')
        dir=s[0]
        fn='c'+s[1]
        
        if ('2004161T131401_2004162T074843' in cubes[i]) or ('2004162T075255_2004162T134240' in cubes[i]): continue
        log=getCamInfo(cubes[i])

        parVis=parseLog(log[1])['Geometry']
        parIR=parseLog(log[2])['Geometry']

        parVis['fileName']=fn
        parVis['dir']=dir
        parIR['fileName']=fn
        parIR['dir']=dir

        #if parIR['SubSpacecraftLongitude']<>None:
        #    print parIR
        #continue

        visObj.append(parVis)
        IRObj.append(parIR)

        
        print >> camInfoLog,log[0]
        for j in range(len(log[1])):
            print >> camInfoLog,log[1][j],
        print >> camInfoLog
        for j in range(len(log[2])):
            print >> camInfoLog,log[2][j],
        print >> camInfoLog,'\n\n'
    camInfoLog.close()
    camInfoPickle=open('/data/VIMS/covims_0004/procdata/camInfo.pick','w+')
    pick.dump([visObj,IRObj],camInfoPickle)
    camInfoPickle.close()
    """
    
    ###parse the camera info output
    #(camInfoVis,camInfoIR)=parseCamInfoLog()
    #camInfoPickle=open('/data/VIMS/covims_0004/procdata/camInfo.pick','w+')
    #pick.dump([camInfoVis,camInfoIR],camInfoPickle)
    #camInfoPickle.close()

    
                    
    with open('/data/VIMS/covims_0004/procdata/camInfo.pick') as camInfoPick:
        [camInfoVis,camInfoIR]=pick.load(camInfoPick)

    o=[]
    for i in range(len(camInfoIR)):
        print camInfoIR[i].keys()
        sys.exit()
        if camInfoIR[i]['SampleResolution']<>None:
            o.append([camInfoIR[i]['SampleResolution'],camInfoIR[i]['dir'],camInfoIR[i]['fileName']])
    o.sort(cmp0)
    for i in range(len(o)):
        print o[i][1],o[i][2],o[i][0]

        
