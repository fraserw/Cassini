import pylab as pyl,numpy as num,pickle as pick,sys

with open('craters/craters_x4_0.pickle') as han:
    (craters,fracCurves,cratMids,taken)=pick.load(han)

#cratersArr=num.array(craters)
#args=num.argsort(cratersArr[:,2])
#cratersArr=cratersArr[args]
#for i in range(len(args)):
#    print cratersArr[i]
#sys.exit()
cratRadss=num.concatenate([num.arange(2.,4.,0.2),num.arange(4.,8.,0.4),num.arange(8.,15.,0.7),num.arange(15.,25.,2.5)])
for i in range(len(craters)):
    print craters[i]

    #suggested algorithm
    #search for an inflection from the minimum and going smaller
    #if the peak is higher than that at 0, bingo
    #otherwise take the 0.
    #then search for a peak at larger radii than the minimum we started from.
    #If the peak is higher than the one we started from, take the peak+1 instead, unless the peak is the maximum point

    arg=num.argmin(fracCurves[i])
    for j in range(arg,0,-1):
        if fracCurves[i][j-1]<fracCurves[i][j]:
            print cratRadss[j+1]
            if fracCurves[i][j]>=fracCurves[i][0]:
                selRad=cratRadss[j+1]
            else:
                selRad=cratRadss[0]
            break

    argmax=num.argmax(fracCurves[i])
    print arg,argmax
    if argmax>arg and argmax<>len(cratRadss)-1:
        selRad=cratRadss[argmax+1]
            
        
    pyl.title(i)
    pyl.plot([craters[i][3],craters[i][3]],[num.min(fracCurves[i]),num.max(fracCurves[i])],'r--')
    pyl.plot([selRad,selRad],[num.min(fracCurves[i]),num.max(fracCurves[i])],'k--')
    pyl.plot(cratRadss,fracCurves[i])
    pyl.show()
