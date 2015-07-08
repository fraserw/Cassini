import pylab as pyl,numpy as num
import pickle as pick

for i in range(8):
    with open('craters/craters_x4_%s.pickle'%(str(i))) as han:
        if i==0:
            [craters,fracCurves,cratMids,taken]=pick.load(han)
            craters=num.array(craters)
        else:
            [c,f,cr,t]=pick.load(han)
            craters=num.concatenate([craters,c])
craters=craters[num.argsort(craters[:,3])]
for i in range(len(craters)):
    print craters[i]
    pyl.scatter((craters[i][1])%360,craters[i][2],s=craters[i][3]*10.)
pyl.show()
