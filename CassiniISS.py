#! /usr/bin/env python

import glob
import commands,os
import sys
import numpy as num, scipy as sci,pickle as pick
from astropy.io import fits as pyf

def findPhoebeCubes(dir):
    cubes=[]    
    files=glob.glob(dir+'/*LBL')
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


def preprocIMGS(cubeName,procDir='/Volumes/data/ISS/coiss_2003/procdata'):
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

    convComm='ciss2isis from=%s.lbl to=%s.cub'%(cube,cube)
    print
    print ' >%s'%(convComm)
    convCommOut=commands.getoutput(convComm)
    print convCommOut
    
    specComm='spiceinit from=%s.cub'%(cube)
    specCommOut=commands.getoutput(specComm)

    print
    print ' >%s'%(specComm)
    print specCommOut

    cisscalComm='cisscal from=%s.cub to=c%s.cub UNITS=I/F'%(cube,cube)
    cisscalOut=commands.getoutput(cisscalComm)

    print
    print ' >%s'%(cisscalComm)
    print cisscalOut

    try:
        os.remove('%s.qub'%cube)
    except:
        pass
    
    os.chdir(cwd)

    log={}
    log[0]=convCommOut
    log[1]=specCommOut
    log[2]=cisscalOut

    return log


def convert2Fits(cubeName,procDir='/Volumes/data/ISS/coiss_2003/procdata'):
    log={}
    cwd=os.getcwd()

    cube=cubeName.split('/')[1].split('.')[0]
    opDir=procDir+'/'+cubeName.split('/')[0]

    os.chdir(opDir)

    convComm='isis2fits from=c%s.cub to=c%s.fits'%(cube,cube)
    print
    print ' >%s '%(convComm)
    out=commands.getoutput(convComm)
    print out
    print
    log[0]=out
    
    os.chdir(cwd)
    return log

if __name__=="__main__":
    #run this from the covims data directory!
    #/Volumes/data/ISS/coiss_2003/data  

    #get Phoebe cubes
    dirs=glob.glob('1465650156_1465674412')
    cubes=[]
    for i in range(len(dirs)):
        out=findPhoebeCubes(dirs[i])
        if out<>None:
            for jj in range(len(out)):
                cubes.append(out[jj])
    
    #pre processing all Phoebe data
    cubeProcLog=open('/Volumes/data/ISS/coiss_2003/procdata/cubeProc.log','w+')    

    for i in range(len(cubes)):
        log=preprocIMGS(cubes[i])
        print >>cubeProcLog,cubes[i]
        for j in [0,1,2]:
            print >>cubeProcLog,log[j]
            print >>cubeProcLog
        print >>cubeProcLog
    cubeProcLog.close()

                
    #convert to fits
    cubeProcLog=open('/Volumes/data/ISSS/coiss_2003/procdata/cubeProc.log','a+')
    print >> cubeProcLog
    for i in range(len(cubes)):
        log=convert2Fits(cubes[i])
        for j in [0,1]:
            print >>cubeProcLog,log[j]
            print >>cubeProcLog
    cubeProcLog.close()
