# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 21:20:41 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/ERA5/o3-1979-2014-JJA-ERA5_2.5.nc')
Ua=xr.open_dataset('D:/bylw/ERA5/uwnd-1979-2014-JJA-ERA5_2.5.nc')
#times=O3.time.data[1979:2014]
levs=O3.level.data[:15]  #UTLS区：250hpa-50hpa
lats=O3.latitude.data[:37] #取北半球纬度
a=O3['o3'].loc[:,:,:0,:]  #选取1979-2014共35年北半球的O3数据
# a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
b=Ua['u'].loc[:,:,:0,:]  #选取1979-2014共35年北半球的Ua数据
# b1=b.loc[b.time.dt.season=='JJA']  #选取夏季数据

#%%分割cell
""" 计算TTu1 """
a10=np.nanmean(a,axis=0)    #对O3求时间平均
a11=np.empty((108,37,37,144))  #对O3求距平分量
for t in range(108):
    a11[t,:,:,:]=a[t,:,:,:]-a10
a20=np.nanmean(a11,axis=3)  #对O3距平分量求纬圈平均
b10=np.nanmean(a,axis=0)    #对Ua求时间平均
b11=np.empty((108,37,37,144))  #对Ua求距平分量
for t in range(108):
    b11[t,:,:,:]=b[t,:,:,:]-b10
b20=np.nanmean(b11,axis=3)  #对Ua距平分量求纬圈平均
ttu1=a20*b20   
 
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
plt.plot([0]*15,levs,linewidth=1,color='k',linestyle='--') 
plt.hlines(50,-65,20,linewidth=1,color='k',linestyle=':')     
plt.plot([0]*15,levs,linewidth=1,color='#EE82EE',label='TFu1',marker='.') 
n.invert_yaxis()          
n.set_yscale('symlog')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
n.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=9.5)
plt.yticks(yticks,ylabel)
plt.xticks(np.arange(-2.5,3,0.5))         
plt.xlim(-2.5,2.5)
n.set_xticklabels(np.arange(-2.5,3,0.5),fontfamily="Times New Roman",fontsize=9.5) 
      
""" 计算TTu2 """      
b21=np.empty((108,37,37,144))  #对Ua距平分量求纬圈距平
for y in range(144):
    b21[:,:,:,y]=b11[:,:,:,y]-b20
ttu2=np.empty((108,37,37,144))  #求TTu2 
for y in range(144):
    ttu2[:,:,:,y]=b21[:,:,:,y]*a20  
ttu20=np.nanmean(ttu2,axis=0)    #对TTu2求时间平均
ttu21=np.empty((37,37,144))   #求纬向散度TTu20 
dlon=2.5*np.pi/180.0   
for x in range(37):  
    dx=6378388.0*np.cos(lats[x]*(np.pi/180.0))*dlon
    for y in range(1,143):
        ttu21[:,x,0]=(ttu20[:,x,1]-ttu20[:,x,0])/dx
        ttu21[:,x,y]=(ttu20[:,x,y+1]-ttu20[:,x,y-1])/(2*dx)     
        ttu21[:,x,143]=(ttu20[:,x,143]-ttu20[:,x,142])/dx
ttu22=ttu21[:,16:24,83:107]
ttu23=np.nanmean(ttu22,axis=(1,2))
# ttu24=xr.DataArray(ttu23,dims=('z', 'x'))
# ttu25=ttu24.rolling(z=3, center=True, min_periods=2).mean()  #对高度求三点滑动平均
# ttu26=np.nanmean(ttu25,axis=1)
plt.plot(-ttu23[:15]*(1e+14),levs,linewidth=1,color='#FFA500',label='TFu2',marker='.')

""" 计算TTu3 """
a21=np.empty((108,37,37,144))  #对O3距平分量求纬圈距平
for y in range(144):
    a21[:,:,:,y]=a11[:,:,:,y]-a20
ttu3=np.empty((108,37,37,144))  #求TTu3 
for y in range(144):
    ttu3[:,:,:,y]=a21[:,:,:,y]*b20  
ttu30=np.nanmean(ttu3,axis=0)    #对TTu2求时间平均
ttu31=np.empty((37,37,144))   #求纬向散度TTu20    
for x in range(37):  
    dx=6378388.0*np.cos(lats[x]*(np.pi/180.0))*dlon
    for y in range(1,143):
        ttu31[:,x,0]=(ttu30[:,x,1]-ttu30[:,x,0])/dx
        ttu31[:,x,y]=(ttu30[:,x,y+1]-ttu30[:,x,y-1])/(2*dx)     
        ttu31[:,x,143]=(ttu30[:,x,143]-ttu30[:,x,142])/dx
ttu32=ttu31[:,16:24,83:107]
ttu33=np.nanmean(ttu32,axis=(1,2))
# ttu34=xr.DataArray(ttu33,dims=('z', 'x'))
# ttu35=ttu34.rolling(z=3, center=True, min_periods=2).mean()  #对高度求三点滑动平均
# ttu36=np.nanmean(ttu35,axis=1)
plt.plot(-ttu33[:15]*(1e+14),levs,linewidth=1,color='#3CB371',label='TFu3',marker='.')         

""" 计算TTu4 """  
ttu4=np.empty((108,37,37,144))  #求TTu4 
for y in range(144):
    ttu4[:,:,:,y]=a21[:,:,:,y]*b21[:,:,:,y]   
ttu40=np.nanmean(ttu4,axis=0)    #对TTu2求时间平均
ttu41=np.empty((37,37,144))   #求纬向散度TTu20    
for x in range(37):  
    dx=6378388.0*np.cos(lats[x]*(np.pi/180.0))*dlon
    for y in range(1,143):
        ttu41[:,x,0]=(ttu40[:,x,1]-ttu40[:,x,0])/dx
        ttu41[:,x,y]=(ttu40[:,x,y+1]-ttu40[:,x,y-1])/(2*dx)     
        ttu41[:,x,143]=(ttu40[:,x,143]-ttu40[:,x,142])/dx
ttu42=ttu41[:,16:24,83:107]
ttu43=np.nanmean(ttu42,axis=(1,2))
# ttu44=xr.DataArray(ttu43,dims=('z', 'x'))
# ttu45=ttu44.rolling(z=3, center=True, min_periods=2).mean()  #对高度求三点滑动平均
# ttu46=np.nanmean(ttu45,axis=1)
plt.plot(-ttu43[:15]*(1e+14),levs,linewidth=1,color='#1E90FF',label='TFu4',marker='.') 
#%%分割cell 
""" 计算TTu总 """
ttu=ttu23+ttu33+ttu43
np.savetxt('tfuls76.csv',ttu,fmt='%.40f',delimiter=',')
plt.plot(-ttu[:15]*(1e+14),levs,linewidth=1,color='k',label='TFu-Total',marker='.')
tu=plt.legend(bbox_to_anchor=(0,-0.25,1,0.2), loc='lower center',mode='expand',ncol=5,prop={'family' : 'Times New Roman', 'size': 9.5})
tu.get_frame().set_linewidth(0.0)
plt.rc('font', family='Times New Roman')
n.set(xlabel=r'浓度变化$\mathrm{-D}$${_x}$${_-}$${_T}$${_F}$$\mathrm{/(10}$${^-}$${^1}$${^4}$ $\mathrm{mol\;·\;s}$${^-}$${^1}$$\mathrm{\;·\;mol}$${^-}$${^1}$$\mathrm{)}$',ylabel=r'气压$\mathrm{/hPa}$',title=r'$\mathrm{(b)}$ 瞬变输送')
plt.show()
print(-ttu[:15])