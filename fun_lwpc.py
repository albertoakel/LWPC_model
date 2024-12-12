#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 13:58:07 2022

@author: akel
"""

import numpy as np
from pysolar.solar import *
import datetime
from mpl_toolkits.basemap import Basemap
import ephem
import math
import subprocess
import matplotlib.pyplot as plt
import datetime
from matplotlib.dates import num2date,datestr2num
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
from astropy.io import fits
from astropy.table import Table



def station():
    tx={'NPM':np.array([21.4,-158]), 
        'NAA': np.array([44.633,-67.283])}
    
    rx={'PLO':np.array([-12.5,-76.8]),
        'CAS':np.array([-31.8,-69.3]),
        'ICA':np.array([-14.0,-75.7]),
        'EACF':np.array([-62.1,-58.4]),
        'TTB':np.array([-1.2,-48.5]),
        'ATI':np.array([-23.18,-46.55])}
    


    return tx,rx

def in_rexp(TX,RX,pts):
    """
    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    pts : integer
        DESCRIPTION.

    Returns
    -------
    d : TYPE
        DESCRIPTION.

    """
    
    TX=TX.upper()
    RX=RX.upper()
    lat,long,d=path(TX,RX,pts) 
    filename='rexp.inp'
    rx_lat=lat[-1]
    rx_lon=long[-1]
    lineTX=["tx-data","     ",TX," ","\n"]
    lineRX=["receivers","     ",str(rx_lat)," ",str(-1.0*rx_lon)," ","\n"]
    dmax=max(d)
    lineRmax=["range-max     ",str(int(dmax))," ","\n"]
    val=round(dmax/pts,1)
    #val=20
    segMax=["lwf-vs-dist 20000"," " ,str(val)," ","\n" ]
    separator =''
    lineTX=separator.join(lineTX)          
    lineRX=separator.join(lineRX)
    lineRmax=separator.join(lineRmax)
    segMax=separator.join(segMax)

    
    arq=open(filename,'w')
    arq.write("case-id     LWPM-REXP \n")
    arq.write("tx          rexp \n")
    arq.write(lineTX)           
    arq.write("ionosphere  range exponential /home/akel/lwpcv21/Profile/test_rexp_py \n")
    arq.write(lineRX)
    arq.write("rx-data     vertical 0 \n") 
    arq.write(lineRmax)                                      
    arq.write("print-swg   3 \n")
    #arq.write("mc mixed \n")
    #arq.write("lwf-vs-dist 20000 25.2 \n")
    arq.write(segMax)
    arq.write("a-noise    ntia  Jul/11/2010 18:01 \n")
    arq.write("start \n") 
    arq.write("quit \n")
    arq.close()
    
    return d

def chi(lat,long,time,**kwargs):
    """
    calculates the zenith angle to array geographic coordinate

    Parameters
    ----------
    lat : geographic latitude ( North is positive, Soult is negative)
    long :geographic longitude(East is positive, west is negative)
        DESCRIPTION.
    time : Date and time in UTC time. example
        time=datetime.datetime(2009,3, 25, 12,0, tzinfo=datetime.timezone.utc)
        DESCRIPTION.

    Returns
    -------
    z : zenith angle

    """
    idx=kwargs.get('idx')

    z=np.zeros(0)
    temp=get_altitude(lat, long, time)-float(90)
    z=np.append(z,temp)
    #print(z)
    
    ##convertendo para valores +_ 
    temp=np.ediff1d(z)
    temp=temp/np.abs(temp)
    #print('-->',len(temp))
    delta=temp[0]
    delta=np.append(delta,temp)
   
    temp2=np.ediff1d(temp)

    indx=np.where(temp2!=0)
    indx=np.asarray(indx)
    delta[np.argwhere(np.isnan(temp))]=1
    if long[0]<long[-1]:
        z=delta*z
    else:
        z=delta*z*(-1.0)

    if idx is not None:
        if indx.size > 0:
            # print('existe 1 valor')
            # print('posição',indx[0,0])
            return z,indx[0,0]
        else:
            #print('Non existe indice')
            indx=0
            return z,indx

    return z


def loadlog_test(filename):
    """
    read infomation logfile
    and return path configuration
    
    Returns
    -------
    V(n,9) :d,lat,long,az,dip,sigma,chi,beta,h'

    """
    #filename='rexp.log'
    skip_line=20
    arq=open(filename,"rt")
    for i in range(skip_line):  
        arq.readline()
    line1=arq.readline()
    temp=np.array(line1.split(),dtype=float)
    v=temp[1:10]
    n=1
    while True:
        n+=1
        line=arq.readline()
        if not line:
            break
        line=line[:-1]
        if line==str(' '):
            #print('Ué')
            break
        temp=np.array(line.split(),dtype=float)
        if len(temp)==9:
            v=np.append(v,temp)
            #print(n,temp[0])
   
    n=int(len(v)/9)
    V=np.reshape(v,(n,9))
    return V

def path(T,R,pt):
    """
    Get Latitude,longitude and distance between Tx-Rx stations
    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    pt : integer
        number of discretized points

    Returns
    -------
    lati : Array float(pt)
        latitude
    longi : Array float(py)
        longitude
    d : Array float(pt)
         distances
    """
    
    mp = Basemap()

    tx,rx=station()
    
    lat=np.array([tx[T][0],rx[R][0]])
    long=np.array([tx[T][1],rx[R][1]])
    
    x1,y1=mp.gcpoints(long[0], lat[0], long[1], lat[1],int(pt))
    
    longi,lati=mp(x1,y1,inverse=True)
    
    R=6372
    lat=np.array(lati)*np.pi/180
    long=np.array(longi)*np.pi/180
    
    d=np.zeros(len(lat)-1)
    for i in range(len(lat)-1):
        k=np.sin(lat[i])*np.sin(lat[int(i+1)])+np.cos(lat[i])*np.cos(lat[i+1])*np.cos(long[i]-long[i+1]) 
        d[i]=R*np.arccos(k)
        
    d=np.insert(d,0,0)
    d=np.cumsum(d)
    return lati,longi,d

def phb(X):
    """
    Define a function for beta and h' for zenith angle
    
    Returns
    -------
    Lat  -latitude
    long -longitude
    d    - distances

    """
    
    #ch=np.array([0.045,1.270,-0.920,3.715,0.295,-1.491,-0.122,0.268])
    ch=np.array([0.045,1.270,-0.774,3.515,0.295,-1.491,-0.122,0.268])
    cb=np.array([0.005,-0.043,0.000,-0.0191,-0.0154])
    H=0.0
    B=0.0
    for i in range(len(ch)):
        H+=ch[i]*X**float(i+1)
    for i in range(len(cb)):
        B+=cb[i]*X**float(i+1)
     
    H=H+70.55# -7
#    B=B+0.395   
    B=B+0.43#+0.18 #- 0.1 #95 #NAA-PLO 25/03/2008

    return H,B

def dchibetah(TX,RX,t,pts):
    """
    

    Parameters
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    t : datetime
        date
    pts :integer
    number of discretized points
.

    Returns
    -------
    h : h' parameter in function
        
    beta : sharpness

    """
    
    lat,long,d=path(TX,RX,pts)
    z,indx=chi(lat,long,t,idx=1)
    
    h=np.zeros(len(z))
    beta=np.zeros(len(z))
    
    for i in range(len(z)):
        if z[i]>-80.0 and z[i]<80:
            X=z[i]*np.pi/180
            h[i],beta[i]=phb(X)
        else:
            h[i]=87
            beta[i]=0.5
            
    if indx!=0:
        c1=indx+1
        c2=indx+2
        dh=np.abs(h[c2]-h[c1])
        db=np.abs(beta[c2]-beta[c1])
        
        if h[c2]>h[c1]:
            h[c2:pts]=h[c2:pts]-dh
        else:
            h[c2:pts]=h[c2:pts]+dh
        
        if beta[c2]>beta[c1]:
            beta[c2:pts]=beta[c2:pts]-db
        else:
            beta[c2:pts]=beta[c2:pts]+db
        
        
    else:
        h=h
    return h,beta

def obsc(time,lat,long):
    """
    

    Parameters
    ----------
    time : TYPE
        DESCRIPTION.
    lat : TYPE
        DESCRIPTION.
    long : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    gatech = ephem.Observer()

    gatech.lon=str(long)
    gatech.lat =str(lat) 

    sun, moon = ephem.Sun(), ephem.Moon()

    results=[]

    gatech.date= (ephem.date(ephem.date(time)))
    sun.compute(gatech)
    moon.compute(gatech)
    r_sun=sun.size/2
    r_moon=moon.size/2
    s=math.degrees(ephem.separation((sun.az, sun.alt), (moon.az, moon.alt)))*60*60

    if s<(r_moon+r_sun):
        lunedelta=0.25*math.sqrt((r_sun+r_moon+s)*(r_moon+s-r_sun)*(s+r_sun-r_moon)*(r_sun+r_moon-s))
    else: 
        lunedelta=None
        percent_eclipse=0
    if lunedelta: 
        lune_area=2*lunedelta + r_sun*r_sun*(math.acos(((r_moon*r_moon)-(r_sun*r_sun)-(s*s))/(2*r_sun*s))) - r_moon*r_moon*(math.acos(((r_moon*r_moon)+(s*s)-(r_sun*r_sun))/(2*r_moon*s)))
        percent_eclipse=(1-(lune_area/(math.pi*r_sun*r_sun)))*100 # Calculate percentage of sun's disc eclipsed using lune area and sun size
    results.append([gatech.date.datetime(),s,sun.size,moon.size,lune_area if lunedelta else 0, percent_eclipse]) 

    return results[0][5]

def run_rexp(TX,RX,t,*k):
    """
    

    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    t : datetime
        date and time
    *k : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pts=80
    in_rexp(TX,RX,pts)
    lat,long,d=path(TX,RX,pts) #obter lat/long
    z=chi(lat,long,t) 
    h,beta=dchibetah(TX,RX,t,pts) 
    time=(t.timetuple()[0],t.timetuple()[1],t.timetuple()[2],t.timetuple()[3],t.timetuple()[4])
    if len(k)!=0:
        OBC=np.zeros(len(lat)) 
        for l in range(len(lat)):
            OBC[l]=obsc(time,lat[l],long[l]) #
        h=h+1*(88-h)*OBC*k*0.01

        
    C=np.stack((d.astype(int), beta,np.round(h,1)),axis=-1)
    np.savetxt('/home/akel/lwpcv21/Profile/test_rexp_py.ndx', C,fmt='%16.2f')
    subprocess.call("./REXP", shell=True)
    out=np.loadtxt('out.dat')
    out=np.delete(out,-1,0)
    return d,h,beta,out,z

def run_rexptime(TX,RX,dat_i,dat_f,dt,**kwargs):
    """
    

    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    dat_i :datetime
        Initial date and time
    dat_f : datetime
        final date and time
    dt : float
        interval in minute
    **kwargs :k
       k=1 inclue eclipse effect

    Returns
    -------
    None.

    """
    K=kwargs.get('k')
    pts=80
    in_rexp(TX,RX,pts) #write inputfile
    lat,long,d=path(TX,RX,pts) #obter lat/long


    if K is not None:
        print("Rum with Eclipse effective")

    delta = datetime.timedelta(minutes=dt)
    t =num2date( drange(dat_i, dat_f, delta)) #rangedate
    amp=np.zeros(0)
    pha=np.zeros(0)
    
    
    for i in range(len(t)):
        z=chi(lat,long,t[i]) 
        h,beta=dchibetah(TX,RX,t[i],pts)
#        h[33:200]=h[33:200]-np.abs(h[32]-h[33])

        time=(t[i].timetuple()[0],t[i].timetuple()[1],t[i].timetuple()[2],t[i].timetuple()[3],t[i].timetuple()[4])
#        OBC=np.zeros(len(lat))
#        for l in range(len(lat)):
#            OBC[l]=obsc(time,lat[l],long[l])
#   K-adjustment factor
#        h=h+1*(88-h)*OBC/100*K   
        C=np.stack((d.astype(int), beta,np.round(h,1)),axis=-1)
        np.savetxt('/home/akel/lwpcv21/Profile/test_rexp_py.ndx', C,fmt='%16.2f') #write geometry h' and beta
        subprocess.call("./REXP", shell=True) #Run 
        out=np.loadtxt('out.dat') #load amp and phase
        amp=np.append(amp,out[-2,1])
        pha=np.append(pha,out[-2,2])
        print(t[i].strftime("%H:%M"),'-->',out[-2,1],out[-2,2])


    return t,amp,pha

def readfit(filename):
    """
    Read fit vlf fit data

    Parameters
    ----------
    filename : STR
        DESCRIPTION.

    Returns
    -------
    C : Array[n,3]
        
    time, amplitude and fase
    

    """
    
    
    temp=fits.open(filename)    
    M=Table(temp[0].data)
    t=M['col0']
    amp=M['col7']
    pha=M['col6']
    C=np.stack((t, amp,pha),axis=-1)
    return C

def mappath(TX,RX,*t):
    """

    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna 
    *t : datetime
        date and time
    Returns
    -------
    None.

    """
    
    
    fig=plt.figure()
    tx,rx=station()   
    lat_d=rx[RX][0]-25
    lat_u=tx[TX][0]+25
    
    lon_l=tx[TX][1]-30
    lon_r=rx[RX][1]+30
   

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
    mp.drawgreatcircle(tx[TX][1],tx[TX][0],rx[RX][1],rx[RX][0],linewidth=2,color='r')
    
    mp.drawparallels(np.arange(-80,90,10),labels=[1,1,0,1])
    mp.drawmeridians(np.arange(-180,180,30),labels=[1,1,0,1])
    title=[TX,'-',RX]
    separator =''
    title=separator.join(title)  
    ax.set_title(title)

    LAT=np.linspace(-90,90,1800)
    LON=np.linspace(-180,180,1800)
    lat, lon = np.meshgrid(LAT, LON)
    mp.fillcontinents()

    if len(t)!=0:     
        Y=t[0].timetuple()[0]
        M=t[0].timetuple()[1]
        D=t[0].timetuple()[2]
        H=t[0].timetuple()[3]
        mn=t[0].timetuple()[4]
        
        t=datetime.datetime(Y,M,D,H,mn,tzinfo=datetime.timezone.utc);

        Z=get_altitude(lat,lon,t)-float(90)
        Z=np.round(Z,2)
        ind90=np.where(Z==-90)
    
        x90,y90 = mp(lon[ind90], lat[ind90])
        mp.plot(x90,y90)
        CS=mp.nightshade(t,alpha=0.1,delta=0.1,color='r',zorder=-800000)
 
    plt.show()
    return
def run_rexp_b(TX,RX,t,*k):
    """
    

    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    t : datetime
        date and time
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    pts=200
    in_rexp(TX,RX,pts)
    lat,long,d=path(TX,RX,pts) #obter lat/long
    z=chi(lat,long,t) 
    h,beta=dchibetah(TX,RX,t,pts) 
    time=(t.timetuple()[0],t.timetuple()[1],t.timetuple()[2],t.timetuple()[3],t.timetuple()[4])
    if len(k)!=0:
        OBC=np.zeros(len(lat)) 
        for l in range(len(lat)):
            OBC[l]=obsc(time,lat[l],long[l]) #
        h=h+1*(88-h)*OBC*k*0.01

        
    C=np.stack((d.astype(int), beta,np.round(h,1)),axis=-1)
    np.savetxt('/home/akel/lwpcv21/Profile/test_rexp_py.ndx', C,fmt='%16.2f')
    subprocess.call("./REXP", shell=True)
    out=np.loadtxt('out.dat')
    out=np.delete(out,-1,0)
    return d,h,beta,out,z

def dchibetah2(TX,RX,t,dh,db,pts):
    """
    

    Parameters
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    t : datetime
        date
    pts :integer
    number of discretized points
.

    Returns
    -------
    h : h' parameter in function
        
    beta : sharpness

    """
    
    lat,long,d=path(TX,RX,pts)
    z,indx=chi(lat,long,t,idx=1)
    
    h=np.zeros(len(z))
    beta=np.zeros(len(z))
    
    for i in range(len(z)):
        if z[i]>-80.0 and z[i]<80:
            X=z[i]*np.pi/180
            h[i],beta[i]=phb2(X,dh,db)
        else:
            h[i]=87
            beta[i]=0.5
            
    if indx!=0:
        c1=indx+1
        c2=indx+2
        dh=np.abs(h[c2]-h[c1])
        db=np.abs(beta[c2]-beta[c1])
        
        if h[c2]>h[c1]:
            h[c2:pts]=h[c2:pts]-dh
        else:
            h[c2:pts]=h[c2:pts]+dh
        
        if beta[c2]>beta[c1]:
            beta[c2:pts]=beta[c2:pts]-db
        else:
            beta[c2:pts]=beta[c2:pts]+db
        
        
    else:
        h=h
    return h,beta

def phb2(X,dh,db):
    """
    Define a function for beta and h' for zenith angle
    
    Returns
    -------
    Lat  -latitude
    long -longitude
    d    - distances

    """
    
    #ch=np.array([0.045,1.270,-0.920,3.715,0.295,-1.491,-0.122,0.268])
    ch=np.array([0.045,1.270,-0.774,3.515,0.295,-1.491,-0.122,0.268])
    cb=np.array([0.005,-0.043,0.000,-0.0191,-0.0154])
    H=0.0
    B=0.0
    for i in range(len(ch)):
        H+=ch[i]*X**float(i+1)
    for i in range(len(cb)):
        B+=cb[i]*X**float(i+1)
     
    H=H+70.55+dh
#    B=B+0.395   
    B=B+0.43+db #- 0.1 #95 #NAA-PLO 25/03/2008

    return H,B
def run_rexp_delta(TX,RX,t,dh,db,*k):
    """
    

    Parameters
    ----------
    TX : STR
        transmitting antenna 
    RX : STR
        receiving antenna
    t : datetime
        date and time
    *k : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    pts=80
    in_rexp(TX,RX,pts)
    lat,long,d=path(TX,RX,pts) #obter lat/long
    z=chi(lat,long,t) 
    h,beta=dchibetah2(TX,RX,t,dh,db,pts) 
    time=(t.timetuple()[0],t.timetuple()[1],t.timetuple()[2],t.timetuple()[3],t.timetuple()[4])
    if len(k)!=0:
        OBC=np.zeros(len(lat)) 
        for l in range(len(lat)):
            OBC[l]=obsc(time,lat[l],long[l]) #
        h=h+1*(88-h)*OBC*k*0.01

        
    C=np.stack((d.astype(int), beta,np.round(h,1)),axis=-1)
    np.savetxt('/home/akel/lwpcv21/Profile/test_rexp_py.ndx', C,fmt='%16.2f')
    subprocess.call("./REXP", shell=True)
    out=np.loadtxt('out.dat')
    out=np.delete(out,-1,0)
    return d,h,beta,out,z