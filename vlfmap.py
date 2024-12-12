
"""
Created on Tue Apr 26 13:42:11 2022

@author: akel
"""

from mpl_toolkits.basemap import Basemap
import numpy as np
import matplotlib.pyplot as plt
import fun_lwpc as lwpc
from datetime import datetime


 
plt.close('all')

def mappath(TX,RX):
    fig=plt.figure()
    tx,rx=lwpc.station()   
    lat_d=rx[RX][0]-25
    lat_u=tx[TX][0]+25
    
    lon_l=tx[TX][1]-30
    lon_r=rx[RX][1]+30
    print(lon_r)

    ax=fig.add_axes([0.1,0.1,0.8,0.8])
    mp = Basemap(urcrnrlat=lat_u,llcrnrlat=lat_d,llcrnrlon=lon_l,urcrnrlon=lon_r,\
                  rsphere=(6378137.00,6356752.3142),\
                  resolution='l',\
                  lat_0=10.,lon_0=-10.,lat_ts=10.)
    #mp=Basemap(projection='mill',lon_0=-120,lat_0=40)
        
        
     
    xtx1,ytx1 = mp(tx[TX][1], tx[TX][0])
    xrx1,yrx1 = mp(rx[RX][1], rx[RX][0])
    mp.plot(xtx1,ytx1,marker='D',color='k')
    mp.plot(xrx1,yrx1,marker='D',color='b')
    
    rx[RX][0]=rx[RX][0]
    print(rx[RX][0],rx[RX][1])
   # mp.drawgreacircle(tx[TX][1],tx[TX][0],rx[RX][1],rx[RX][0],linewidth=2,color='r')
    mp.drawgreatcircle(tx[TX][1],tx[TX][0],rx[RX][1],rx[RX][0],linewidth=2,color='r')
    
    mp.fillcontinents()
    mp.drawparallels(np.arange(-80,90,10),labels=[1,1,0,1])
    mp.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
    title=[TX,'-',RX]
    separator =''
    title=separator.join(title)  
    ax.set_title(title)
    date = datetime(2008,6, 25,22,0);

#    CS=mp.nightshade(date,alpha=0.2,color='m')
    
    
    plt.show()
    return



mappath('NPM','NAA')



