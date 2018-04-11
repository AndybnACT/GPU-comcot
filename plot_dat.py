from mpl_toolkits.basemap import Basemap
from sys import argv
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as npy
import struct
from os import listdir
class COMCOTdat(object):
    """COMCOTdat
    Import comcot data from .dat files

    """
    def __init__(self, BASE_DIR, Layer=1):
        super(COMCOTdat, self).__init__()
        if BASE_DIR[-1] != '/':
            BASE_DIR+='/'
        # import data
        self._x = self.ReadFile(BASE_DIR+ "layer%02d_x.dat" % Layer)
        self._y = self.ReadFile(BASE_DIR+ "layer%02d_y.dat" % Layer)
        self._dim = [len(self._y), len(self._x)]
        self._z = self.ReadFile(BASE_DIR+ "layer%02d.dat" % Layer, dim = True)
        #self._zmax = self.ReadFile(BASE_DIR+ "zmax_layer%02d.dat" % Layer, dim = True)
        self._etaname = [i for i in listdir(BASE_DIR) if i[:4]=='z_%02d' % Layer]
        self._etaname.sort()
        self._eta = {}
        for timestrig in self._etaname:
            minute = int(timestrig[5:-4])
            masked = npy.ma.masked_where(self._z < 0, self.ReadFile(timestrig, dim=True, binary=True))
            self._eta.update({minute:masked})
        # self._arri = self.ReadFile(BASE_DIR+ "arrival_layer%02d.dat" % Layer,dim = True)


        [self._xx, self._yy] = npy.meshgrid(self._x,self._y)

    def ReadFile(self,filename,dim=False , binary=False):
        print("Reading File %s" % filename)
        try: # read from file
            if binary:
                with open(filename,'rb') as FILE:
                    size = self.dim[0]*self.dim[1]
                    data = struct.unpack('f'*size, FILE.read(4*size))
                    if dim:
                        return npy.reshape(data,self.dim)
                    else:
                        return data
            else:
                with open(filename,'r') as FILE:
                    text = FILE.read()
                    #parse to list
                    data = [float(value) for value in text.split()]
                    if dim:
                        return npy.reshape(data,self.dim)
                    else:
                        return data
        except FileNotFoundError as e:
            print("Error, File < %s > Not Found" % e.filename)
            raise


    @property
    def dim(self):
        return self._dim
    @property
    def XGRID(self):
        return self._xx
    @property
    def YGRID(self):
        return self._yy

    @property
    def ZMAX(self):
        return self._zmax
    @property
    def ARRIVAL(self):
        return self._arri
    @property
    def ETA(self):
        return self._eta

class COMCOTGlobalmap(COMCOTdat):
    def __init__(self, BASE_DIR, Layer=1):
        super(COMCOTGlobalmap, self).__init__(BASE_DIR, Layer=1)
        self._cmax = 0.5#npy.max(self.ETA[0])/2
        for time in self.ETA:
            self.Eta_Map(time,2)
            plt.savefig("plot_z_01_%05d.png" %time)
            plt.close()
    def Eta_Map(self, step, timestep):
        print('processing %dth step plot' % step)
        m = self.Global_MAP()
        pc = m.pcolor(self.XGRID, self.YGRID, self.ETA[step], cmap='RdBu_r', vmin=-self._cmax, vmax=self._cmax)
        m.colorbar(pc, location="right")
        step = float(step)/(60.0/timestep)
        plt.title("Tsunami Wave Propagation at Time = %5.3f min" % step)

    def Global_MAP(self):
        fig = plt.figure(figsize=(18,12))
        fig.add_axes([0.05,0.05,0.9,0.9])
        ax = plt.gca()
        m = Basemap(projection='cyl',llcrnrlat=0.0,urcrnrlat=self.YGRID[-1,0],\
                llcrnrlon=self.XGRID[0,0],urcrnrlon=255,resolution='i', ellps="WGS84", ax=ax)
        # h = 5000.
        # m = Basemap(projection='nsper',lon_0=139,lat_0=35,
        #             satellite_height=h*1000.,resolution='l', ax=ax)
        m.drawcoastlines()
        m.fillcontinents(color='coral', lake_color='aqua')
        m.drawmapboundary(fill_color='aqua')
        # draw parallels and meridians.
        m.drawparallels(npy.arange(-90.,91.,20.),labels=[1,0,0,0])
        m.drawmeridians(npy.arange(-180.,181.,15.),labels=[0,0,0,1])
        return m
    
if __name__ == '__main__':
    test = COMCOTGlobalmap(argv[1])
