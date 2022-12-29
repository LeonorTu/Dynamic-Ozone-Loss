# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 21:57:40 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/ERA5/o3-1979-2014-JJA-ERA5_2.5.nc')
Va=xr.open_dataset('D:/bylw/ERA5/vwnd-1979-2014-JJA-ERA5_2.5.nc')
#times=O3.time.data[1979:2014]
levs=O3.level.data[:15]  #UTLS区：250hpa-50hpa
lats=O3.latitude.data[36:] #取北半球纬度
stu=np.loadtxt('sfuls76.csv',delimiter=',')
ttu=np.loadtxt('tfuls76.csv',delimiter=',')
stv=np.loadtxt('sfvls76.csv',delimiter=',')
ttv=np.loadtxt('tfvls76.csv',delimiter=',')
sttt=stu-stv
tttt=ttu-ttv
utt=stu+ttu
vtt=-stv-ttv
total=stu-stv+ttu-ttv
config = {
    "font.family":'serif',
    "font.size": 10.5,
    "mathtext.fontset":'stix',
    "font.serif": ['STZhongsong'],
}
plt.rcParams.update(config) 
fig=plt.figure(1,dpi=800)   #画图
plt.rcParams['figure.figsize'] = (4,3) # 单位是inches
n=fig.add_axes([0, 0, 1, 1])    
plt.hlines(50,-620,460,linewidth=1,color='k',linestyle=':')       
plt.plot([0]*15,levs,linewidth=1,color='k',linestyle='--') 
plt.plot(-sttt[:15]*(1e+14),levs,linewidth=1,color='#EE82EE',label='SF-Total',marker='.') 
plt.plot(-tttt[:15]*(1e+14),levs,linewidth=1,color='#FFA500',label='TF-Total',marker='.') 
plt.plot(-utt[:15]*(1e+14),levs,linewidth=1,color='#3CB371',label='u-Total',marker='.') 
plt.plot(-vtt[:15]*(1e+14),levs,linewidth=1,color='#1E90FF',label='v-Total',marker='.') 
plt.plot(-total[:15]*(1e+14),levs,linewidth=1,color='k',label='Total',marker='.') 
n.invert_yaxis()          
n.set_yscale('log')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
n.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=9.5)
plt.yticks(yticks,ylabel)
plt.xticks(np.arange(-140,81,20))         
plt.xlim(-140,80)
n.set_xticklabels(np.arange(-140,81,20),fontfamily="Times New Roman",fontsize=9.5)   
tu=plt.legend(bbox_to_anchor=(0,-0.25,1,0.2), loc='lower center',mode='expand',ncol=5,prop={'family' : 'Times New Roman', 'size': 9.5})
tu.get_frame().set_linewidth(0.0)
n.set(xlabel=r'浓度变化$\mathrm{-D}$$\mathrm{/(10}$${^-}$${^1}$${^3}$ $\mathrm{mol\;·\;s}$${^-}$${^1}$$\mathrm{\;·\;mol}$${^-}$${^1}$$\mathrm{)}$',ylabel=r'气压$\mathrm{/hPa}$',)
plt.show()
# print(utt)
# print(vtt)
print(-sttt[:15])