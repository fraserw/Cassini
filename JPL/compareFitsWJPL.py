import pickle as pick
import numpy as np,pylab as pyl
import jpleph

jpl = eph.eph()


with open('fitsTable.dat') as han:
    data = han.readlines()

imageNames, sampleResolutionnot, long_o_not,lat_o_not,az_o_not,long_s,lat_s,distancenot,offXV,offYV,offXI,offYI,offsetVel,offsetVelAngle = [],[],[],[],[],[],[],[],[],[],[],[],[],[]
utc = []
for i in range(len(data)):
    s = data[i].split()
    imageNames.append(s[0])
    sampleResolutionnot.append(float(s[1]))
    long_o_not.append(float(s[2]))
    lat_o_not.append(float(s[3]))
    az_o_not.append(float(s[4]))
    long_s.append(float(s[5]))
    lat_s.append(float(s[6]))
    distancenot.append(float(s[7]))
    offXV.append(float(s[8]))
    offYV.append(float(s[9]))
    offXI.append(float(s[10]))
    offYI.append(float(s[11]))
    offsetVel.append(float(s[12]))
    offsetVelAngle.append(float(s[13]))

distancenot = np.array(distancenot)
long_s = np.array(long_s)
long_o = np.array(long_o_not)
lat_s = np.array(lat_s)
lat_o = np.array(lat_o_not)


jplDists = []
JDs = []
jplLong_ss = []
jplLat_ss = []
jplLong_os = []
jplLat_os = []
for i in range(len(imageNames)):
    with open('/data/VIMS/covims_0004/procdata/%s_vis.campt.pickle'%(imageNames[i])) as han:
        latLongObjVis=pick.load(han)
    junk = []
    for j in range(len(latLongObjVis)):
        for k in range(len(latLongObjVis[j])):
            if latLongObjVis[j][k]['SpacecraftPosition']<>None:
                junk.append(latLongObjVis[j][k]['UTC'])

    (JD,jplDist,jplLong_o,jplLat_o,jplLong_s,jplLat_s) = jpl(junk[0])
    JDs.append(JD)
    jplDists.append(jplDist)
    jplLong_ss.append(jplLong_s)
    jplLat_ss.append(jplLat_s)
    jplLong_os.append(jplLong_o)
    jplLat_os.append(jplLat_o)
    print imageNames[i],distancenot[i],jplDists[i]


jplLong_ss = np.array(jplLong_ss)
jplLong_os = np.array(jplLong_os)
jplLat_ss = np.array(jplLat_ss)
jplLat_os = np.array(jplLat_os)
jplDists = np.array(jplDists)

diff = (long_o-jplLong_os)
print diff
print np.mean(diff),np.std(diff)
pyl.plot(JDs,diff)

pyl.show()

#offsets
long_o = 1.609
long_s = 3.578
lat_s = 7.911
fractionalDistance = 1.0499
