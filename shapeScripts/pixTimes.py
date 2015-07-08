import numpy as num
from lineFitUtils import poly
import sys


###this code appears to work fine for all vims images, even in 1x1 binned visual channel
def getPixTimes(Cass,A,B,Tnot,inot,jnot):
    #VIMS scans along samples from 1 to 35, then increases the line
    VIMSexps=[]
    for i in range(len(Cass)):
        if Cass[i][5]==Cass[i-1][5] and Cass[i][6]==Cass[i-1][6]+1:
            VIMSexps.append(Cass[i][3]-Cass[i-1][3])

    #time to take a pixel in the same line
    singlePixExpTime=num.median(VIMSexps)


    #select the most common sample pixel #
    w=num.where(Cass[:,6]== num.bincount(Cass[:,6].astype('int')).argmax())
    x=Cass[w][:,3]
    
    #time to swap from one line to the next
    lineSwitchTime=num.median((x[1:]-x[:len(x)-1]))-B*singlePixExpTime

    #setup the time array. Pick as Tnot the first time in the Cass array
    pixelTimes=num.zeros([A,B]).astype(num.float64)
    
    t=0.0
    for i in range(A):
        for j in range(B):
            pixelTimes[i,j]=t
            t+=singlePixExpTime
        t+=lineSwitchTime
    
    pixelTimes+=Tnot-pixelTimes[inot,jnot]

    return pixelTimes


def getVels(Cass,pixelTimes,Xnot,Ynot,Znot,inot,jnot,polyOrder=2):#velocity polynomial fitting order

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
    return (spaceCraftVectornot,deltaX,deltaY,deltaZ)
