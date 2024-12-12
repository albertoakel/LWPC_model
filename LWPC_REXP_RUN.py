
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:28:09 2022
run to LWPC in Ranbge exponetial mode
define 
@author: akel
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from numpy import loadtxt
import subprocess
from pysolar.solar import *
from matplotlib.dates import (YEARLY, DateFormatter,
                              rrulewrapper, RRuleLocator, drange)
from matplotlib.dates import num2date,datestr2num
import datetime
import fun_lwpc as lwpc
from pandas.plotting import andrews_curves






plt.close('all')

#V=lwpc.loadlog_test('lwpm.log') #NAA-PLO 2009/03/25 0130
# d0=V[:,0]
# lat=V[:,1]
# long=V[:,2]*-1
#chi=V[:,6]
# b=V[:,7]
# h=V[:,8]

dat1 = datetime.datetime(2008, 3, 25,12,0,tzinfo=datetime.timezone.utc);
dat2 = datetime.datetime(2008, 3, 25,22,0,tzinfo=datetime.timezone.utc);
delta = datetime.timedelta(minutes=60)
t =num2date( drange(dat1, dat2, delta))



#z=np.zeros((len(t),600))

#t = datetime.datetime(2020,12,14,17,10, tzinfo=datetime.timezone.utc)
amp=np.zeros(0)
pha=np.zeros(0)

#z=np.zeros(0)
pts=200

for i in range(len(t)):
    lat,long,d=lwpc.path('NAA','PLO',pts) #get  lat/long
    z=lwpc.chi(lat,long,t[i])
    
    h,beta=lwpc.dchibetah('NAA','PLO',t[i],pts) #evaluation chi e beta
    time=(t[i].timetuple()[0],t[i].timetuple()[1],t[i].timetuple()[2],t[i].timetuple()[3],t[i].timetuple()[4])
    obsc=np.zeros(len(lat))
    # for l in range(len(lat)):
    #     obsc[l]=lwpc.obsc(time,lat[l],long[l])
    
    # k=0.0
    # h=h+1*(88-h)*obsc/100*k   
    C=np.stack((d, beta,h),axis=-1)
    np.savetxt('/home/akel/lwpcv21/Profile/test_rexp_py.ndx', C,fmt='%-7.2f')
    subprocess.call("./REXP", shell=True)
    out=np.loadtxt('out.dat')
    d=out[:,0]
    #temp=out[-1,1]
    print(t[i].strftime("%H:%M"),'-->',out[-2,1],out[-2,2])
    #print(t[i].timetuple()[3],t[i].timetuple()[4],'--->',out[-1,1])
    amp=np.append(amp,out[-2,1])
    pha=np.append(pha,out[-2,2])




        


dat_lim = datetime.datetime(2008, 3, 25,21,30,tzinfo=datetime.timezone.utc);
minor_ticks =drange(dat1, dat2, datetime.timedelta(minutes=5)) #


#t=np.linspace(12,21.75,len(amp))
fig1, ax1 = plt.subplots()
     
ax1.plot(t,amp)
ax1.set_xticks(minor_ticks, minor=True)

formatter = DateFormatter('%H:%M')
ax1.xaxis.set_major_formatter(formatter)
plt.xlim(dat1,dat_lim)

fig2, ax2 = plt.subplots()        
ax2.plot(t,pha)
ax2.set_xticks(minor_ticks, minor=True)
formatter = DateFormatter('%H:%M')
ax2.xaxis.set_major_formatter(formatter)
plt.xlim(dat1,dat_lim)
#plt.ylim(-100,0)



# k=0.1
# h1=h+0*(88-h)*obsc/100*k  #eclipse effect
#fig1, ax3 = plt.subplots()
# ax3.plot(d,h,'b',d,h1,'--')
# ax3.set_ylabel('h+obsc')
# ax3.set_xlabel('Distancia') 
#plt.ylim([70,78])


# plt.xlim([0,7000])    
# plt.ylim([0.34,0.5])
    
    

    
    
    







