import numpy as np

class eph:
    def __init__(self):
        with open('/Users/fraserw/git/Cassini/JPL/jpl_eph.dat') as han:
            data = han.readlines()

        self.CP_dist = []
        self.utc_strings = []
        self.JD = []
        self.long_s = []
        self.lat_s = []
        self.long_o = []
        self.lat_o = []
        for i in range(len(data)):
            s = data[i].split()
            self.utc_strings.append(s[0].replace('-Jun-','-06-')+'T'+s[1])
            self.CP_dist.append(float(s[7]))
            self.JD.append(float(s[2]))
            self.long_o.append(float(s[3]))
            self.lat_o.append(float(s[4]))
            self.long_s.append(float(s[5]))
            self.lat_s.append(float(s[6]))
        self.CP_dist = np.array(self.CP_dist) * 149597870700.0/1.e3
        self.utc_strings = np.array(self.utc_strings)
        self.JD = np.array(self.JD)
        self.long_s = 360.0-np.array(self.long_s)
        self.lat_s = np.array(self.lat_s)
        self.long_o = 360.0-np.array(self.long_o)
        self.lat_o = np.array(self.lat_o)


        self.long_o_corr = 1.609
        self.long_s_corr = 3.578
        self.lat_s_corr = 7.911
        self.fractionalDistance_corr = 1.0499

    def __call__(self,time,includeCorrection = True):
        t = time[:16]
        w = np.where(self.utc_strings == t)
        if includeCorrection:
            return (self.JD[w][0],
                    self.CP_dist[w][0]*self.fractionalDistance_corr,
                    self.long_o[w][0]+self.long_o_corr,
                    self.lat_o[w][0],
                    self.long_s[w][0]+self.long_s_corr,
                    self.lat_s[w][0]+self.lat_s_corr)
        else:
            return (self.JD[w][0],self.CP_dist[w][0],self.long_o[w][0],self.lat_o[w][0],self.long_s[w][0],self.lat_s[w][0])
