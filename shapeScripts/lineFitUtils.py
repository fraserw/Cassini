#! /usr/bin/env python

import numpy as num
from numpy import linalg
import sys

class poly:
    def __init__(self,x,y,order=2,renormX=False):
        """Create a polynomial object of the form
              Y=sum_i^order c_i*X**i
           using least squares to determine the coefficients c_i

           Fit is done on initialization. Self call of the object
           will return the best-fit Y values.

           Order will determine the order of the polynomial. Default
           is 2, ie., a line.

           if renormX is true, the median of x will be
           removed from the array to prevent numerical
           errors.
        """
        self.renorm=renormX
        if self.renorm:
            self.medx=num.median(x)
            self.x=x-self.medx
        else:
            self.x=x
            self.medx=0.
        self.y=y
        self.order=order

        self._fit()

    def changeOrder(self,newOrder):
        """Reevaluate the least squares fit using the new
           order.
        """
        self.order=newOrder
        self._fit()

    def _fit(self):
        self.A=num.ones([self.x.shape[0],self.order]).astype(self.y.dtype)
        for ii in range(1,self.order):
            self.A[:,ii]=self.x**ii
        At=self.A.transpose()
        Abrack=num.dot(At,self.A)
        AbrackInv=linalg.inv(Abrack)
        self.poly=num.dot(num.dot(linalg.inv(Abrack),At),self.y)

    def _eval(self):
        self.Y=num.dot(self.A,self.poly)

    def __call__(self):
        self._eval()
        return self.Y

    def polyRMS(self):
        """Returns the rms of the LSQ fit
        """
        self._eval()
        return num.sum( (self.Y-self.y)**2)**0.5

    def eval(self,X):
        """Evaluate the polynomial at the provided X values
        """
        A=num.ones([X.shape[0],self.order]).astype(self.y.dtype)
        for ii in range(1,self.order):
            A[:,ii]=(X-self.medx)**ii
        return num.dot(A,self.poly)
    
if __name__=="__main__":
    import scipy as sci
    x=num.linspace(0,100,200)
    y=0.4*x+200.
    y+=sci.randn(len(x))

    obj=poly(x,y,order=2,renormX=False)
    print obj.polyRMS()
    obj.changeOrder(4)
    print obj.polyRMS()
