import matplotlib
matplotlib.use('TkAgg')
from matplotlib import image as mpimg
import pylab as pyl
import numpy as num
from matplotlib.patches import Ellipse

class ev:

    def __init__(self,figure,sp,image):
        self.patchList=[]
        #self.patchRads=[]
        self.sp=sp
        self.image=image

        self.shift_is_held=False
        self.ctrl_is_held=False
        self.press=None

        cidbutton = figure.canvas.mpl_connect('button_press_event', self.onclick)
        cidpick = figure.canvas.mpl_connect('pick_event', self.onpick)
        cidpress = figure.canvas.mpl_connect('key_press_event',self.on_key_press)
        cidrelease = figure.canvas.mpl_connect('key_release_event',self.on_key_release)
        cidmotion = figure.canvas.mpl_connect('motion_notify_event',self.onmotion)

    def loadCraters(self,cratFile='manual.craters'):
        with open(cratFile) as han:
            data=han.readlines()
        for ii in range(len(data)):
            s=data[ii].split()
            (x,y,h)=(float(s[0]),float(s[1]),float(s[2]))
            w=self.getWidthAdjust(y)*h
            ellipse=Ellipse((x,y), w,h, angle=0.,
            facecolor="blue", edgecolor="red",zorder=10, linewidth=2, alpha=0.55,
            picker=True)
            self.patchList.append(ellipse)
            self.sp.add_patch(ellipse)
            pyl.draw()
    def getWidthAdjust(self,y):
        a=self.image.shape[0]
        lat=num.arcsin(2*(0.5*a-y)/float(a))
        r=1./num.cos(lat)
        return r

    def onclick(self,event):

        posx=event.xdata
        posy=event.ydata
        if posx<>None and posy<>None and event.button==1 and not self.shift_is_held and not self.ctrl_is_held:

            cratRadius=50.
            r=cratRadius*self.getWidthAdjust(posy)


            ellipse=Ellipse((posx,posy), r,cratRadius, angle=0.,
            facecolor="blue", edgecolor="red",zorder=10, linewidth=2, alpha=0.55,
            picker=True)
            self.sp.add_patch(ellipse)
            self.patchList.append(ellipse)
            #self.patchRads.append(cratRadius)

            pyl.draw()

    def on_key_press(self,event):
        if event.key=='shift':
            self.shift_is_held=True
        if event.key=='control':
            self.ctrl_is_held=True
        if event.key in ['w','W']:
            with open('manual.craters','w+') as han:
                for ii in range(len(self.patchList)):
                    print >>han,self.patchList[ii].center[0],self.patchList[ii].center[1],self.patchList[ii].height

    def on_key_release(self,event):
        if event.key=='shift':
            self.shift_is_held=False
            self.press=None
        if event.key=='control':
            self.ctrl_is_held=False
            self.press=None

    def onpick(self,event):
        if event.mouseevent.button==3and not self.shift_is_held:
            for ii in range(len(self.patchList)):
                if event.artist==self.patchList[ii]:
                    k=ii
                    break
            self.patchList[k].remove()
            pyl.draw()
        if event.mouseevent.button==1 and self.shift_is_held:
            for ii in range(len(self.patchList)):
                if event.artist==self.patchList[ii]:
                    k=ii
                    break

            self.press=(self.patchList[k].center[0],self.patchList[k].center[1],event.mouseevent.xdata,event.mouseevent.ydata,k,self.patchList[k].height,1)

        elif event.mouseevent.button==1 and self.ctrl_is_held:
            for ii in range(len(self.patchList)):
                if event.artist==self.patchList[ii]:
                    k=ii
                    break

            self.press=(self.patchList[k].center[0],self.patchList[k].center[1],event.mouseevent.xdata,event.mouseevent.ydata,k,self.patchList[k].height,3)

    def onmotion(self,event):
        if self.press==None: return
        (x0,y0,xpress,ypress,k,h,button)=self.press

        if button==1:
            dx=event.xdata-xpress
            dy=event.ydata-ypress
            self.patchList[k].center=(x0+dx,y0+dy)
            self.patchList[k].width=h*self.getWidthAdjust(event.ydata)
            #self.patchList[k].width=self.patchRads[k]*self.getWidthAdjust(event.ydata)
            self.patchList[k].figure.canvas.draw()
        elif button==3:
            dx=event.xdata-xpress
            self.patchList[k].height=h*(1.+dx/self.image.shape[1]*5.)
            self.patchList[k].width=h*(1.+dx/self.image.shape[1]*5.)*self.getWidthAdjust(y0)
            self.patchList[k].figure.canvas.draw()


fig=pyl.figure('Phoebe',figsize=(30,15))
sp1=fig.add_subplot(111)

image=mpimg.imread('/data/PhoebePDSMaps/PhoebeFull.png')
print image.shape
im=pyl.imshow(image)

ob=ev(fig,sp1,image)
ob.loadCraters()
pyl.show()