#! /usr/bin/env python

import pickle as pick,numpy as num

def getFit(fn):
    with open(fn) as han:
        (samps,probs)=pick.load(han)
    
    (Y,X,nDim)=samps.shape
    goodSamps=[]
    for i in range(Y):
        for j in range(int(X*0.5),X):
            xx=[]
            for k in range(nDim):
                xx.append(samps[i,j][k])
            xx.append(probs[i,j])
            goodSamps.append(xx[:])
    goodSamps=num.array(goodSamps)
    (l,ll)=goodSamps.shape
    args=num.argsort(goodSamps[:,ll-1])
    goodSamps=goodSamps[args]

    return (goodSamps[l-1],goodSamps)
    
if __name__=="__main__":
    (bp,gs)=getFit('cv1465671822_1.noVel_fit_pickle')
    print bp
    print gs.shape
