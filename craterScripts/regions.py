import numpy as num
import sys

class region:
    def __init__(self,borders,step=0.01):

        self.borders=borders

        minx=num.min(borders[:,0])
        miny=num.min(borders[:,1])
        maxx=num.max(borders[:,0])
        maxy=num.max(borders[:,1])

        x=num.arange(minx,maxx+step,step)
        y=num.zeros((len(x),2))
        y[:,0]+=1000.
        y[:,1]-=1000.

        for ii in range(len(borders)-1):
            m=(borders[ii+1][1]-borders[ii][1])/(borders[ii+1][0]-borders[ii][0])
            b=borders[ii][1]-m*borders[ii][0]
            Y=m*x+b
            w=num.where((Y>=min(borders[ii][1],borders[ii+1][1]))&(Y<=max(borders[ii][1],borders[ii+1][1])))
            for j in w[0]:
                y[j][0]=min(y[j,0],Y[j])
                y[j][1]=max(y[j,1],Y[j])
        if y[0][0]==1000.:
            y=y[1:]
        if y[-1,0]==1000.:
            y=y[:-1]
        self.y=y
        self.x=x

    def __call__(self,x,y):
        if x<self.x[0] or x>self.x[-1]: return False
        arg=num.argmin(num.abs(self.x-x))
        if y>=self.y[arg][0] and y<=self.y[arg][1]:
            return True
        else:
            return False

if __name__=="__main__":
    poorRegion1=region(num.array([[119,-40.7],
    [83.7,-31.6],
    [67.0,-13.8],
    [58.4,-2.5],
    [65.5,18.3],
    [74.9,27.9],
    [90.2,21.8],
    [94.0,7.0],
    [119.0,-40.7]]))
    print poorRegion1(30.,30.)
    print poorRegion1(70.,10.)
