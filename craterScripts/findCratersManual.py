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

        self.cratRad=5.

        cidbutton = figure.canvas.mpl_connect('button_press_event', self.onclick)
        cidpick = figure.canvas.mpl_connect('pick_event', self.onpick)
        cidpress = figure.canvas.mpl_connect('key_press_event',self.on_key_press)
        cidrelease = figure.canvas.mpl_connect('key_release_event',self.on_key_release)
        cidmotion = figure.canvas.mpl_connect('motion_notify_event',self.onmotion)

    def remPatch(self,l):
        if l==0:
            self.patchList[l].remove()
            self.patchList=self.patchList[1:]
        elif l==len(self.patchList)-1:
            self.patchList[l].remove()
            self.patchList=self.patchList[:l]
        else:
            self.patchList[l-1].remove()
            self.patchList=self.patchList[:l-1]+self.patchList[l+1:]

    def loadCraters(self,cratFile='manual.craters'):
        with open(cratFile) as han:
            data=han.readlines()
        crats=[]
        for ii in range(len(data)):
            s=data[ii].split()
            (x,y,h)=(float(s[0]),float(s[1]),float(s[2]))
            crats.append([x,y,h])
        crats=num.array(crats)
        arg=num.argsort(crats[:,2])
        crats=crats[arg]
        for ii in range(len(crats)):
            (x,y,h)=crats[ii]
            w=self.getWidthAdjust(y)*h
            ellipse=Ellipse((x,y), w,h, angle=0.,
            facecolor="blue", edgecolor="red",zorder=len(crats)-ii, linewidth=2, alpha=0.55,
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
        print pyl.get_current_fig_manager().toolbar._active,posx
        if pyl.get_current_fig_manager().toolbar._active==None and posx<>None and posy<>None and event.button==1 and not self.shift_is_held and not self.ctrl_is_held:

            r=self.cratRad*self.getWidthAdjust(posy)


            ellipse=Ellipse((posx,posy), r,self.cratRad, angle=0.,
            facecolor="blue", edgecolor="red",zorder=10, linewidth=2, alpha=0.55,
            picker=True)
            self.sp.add_patch(ellipse)
            self.patchList.append(ellipse)
            #self.patchRads.append(cratRadius)

            pyl.draw()

    def on_key_press(self,event):
        if event.key=='?':
            help_string="""Press and hold shift, then select and drag to move a circle.
            Press and hold control, then move mouse left or right to resize a circle.
            Press 1 to 5 to select a default circle size (1 smallest, 5 largest).
            Press capital D to delete the last created circle.

            Left click to create a new circle.
            Right click on a circle to delete it.

            Press w or W to save the marked crater regions to the file manual.craters
            """
            print help_string
        if event.key=='shift':
            self.shift_is_held=True
        if event.key=='control':
            self.ctrl_is_held=True
        if event.key in ['w','W']:
            with open('manual.craters','w+') as han:
                for ii in range(len(self.patchList)):
                    print >>han,self.patchList[ii].center[0],self.patchList[ii].center[1],self.patchList[ii].height
        if event.key=='D':
            print 'Deleting last created crater.'
            self.remPatch(len(self.patchList))
            #l=len(self.patchList)
            #self.patchList[l-1].remove()
            #self.patchList=self.patchList[:l-1]
            pyl.draw()

        if event.key=='1':
            self.cratRad=1.
            print self.cratRad
        if event.key=='2':
            self.cratRad=6.
            print self.cratRad
        if event.key=='3':
            self.cratRad=10.
            print self.cratRad
        if event.key=='4':
            self.cratRad=30.
            print self.cratRad
        if event.key=='5':
            self.cratRad=50.
            print self.cratRad
        if event.key=='6':
            self.cratRad=70.
            print self.cratRad


    def on_key_release(self,event):
        if event.key=='shift':
            self.shift_is_held=False
            self.press=None
        if event.key=='control':
            self.ctrl_is_held=False
            self.press=None

    def onpick(self,event):
        if event.mouseevent.button==3 and not self.shift_is_held and pyl.get_current_fig_manager().toolbar._active==None:
            for ii in range(len(self.patchList)):
                if event.artist==self.patchList[ii]:
                    k=ii
                    break
            self.remPatch(k)
            #self.patchList[k].remove()
            pyl.draw()

        elif event.mouseevent.button==1 and self.shift_is_held  and pyl.get_current_fig_manager().toolbar._active==None:
            for ii in range(len(self.patchList)):
                if event.artist==self.patchList[ii]:
                    k=ii
                    break

            self.press=(self.patchList[k].center[0],self.patchList[k].center[1],event.mouseevent.xdata,event.mouseevent.ydata,k,self.patchList[k].height,1)

        elif event.mouseevent.button==1 and self.ctrl_is_held and pyl.get_current_fig_manager().toolbar._active==None:
            for ii in range(len(self.patchList)):
                if event.artist==self.patchList[ii]:
                    k=ii
                    break

            self.press=(self.patchList[k].center[0],self.patchList[k].center[1],event.mouseevent.xdata,event.mouseevent.ydata,k,self.patchList[k].height,3)

    def onmotion(self,event):

        if pyl.get_current_fig_manager().toolbar._active==None:
            if self.press==None: return
            (x0,y0,xpress,ypress,k,h,button)=self.press

            if button==1 and pyl.get_current_fig_manager().toolbar._active==None:
                dx=event.xdata-xpress
                dy=event.ydata-ypress
                self.patchList[k].center=(x0+dx,y0+dy)
                self.patchList[k].width=h*self.getWidthAdjust(event.ydata)
                #self.patchList[k].width=self.patchRads[k]*self.getWidthAdjust(event.ydata)
                self.patchList[k].figure.canvas.draw()
            elif button==3 and pyl.get_current_fig_manager().toolbar._active==None:
                dx=event.xdata-xpress
                multi=max((1.+dx/self.image.shape[1]*5.),1.e-3)
                self.patchList[k].height=h*multi
                self.patchList[k].width=h*multi*self.getWidthAdjust(y0)
                self.patchList[k].figure.canvas.draw()

if __name__=="__main__":
    fig=pyl.figure('Phoebe',figsize=(16,8))
    sp1=fig.add_subplot(111)

    image=mpimg.imread('/data/PhoebePDSMaps/PhoebeFull.png')
    print image.shape
    im=pyl.imshow(image,zorder=0)

    ob=ev(fig,sp1,image)
    ob.loadCraters()
    pyl.show()