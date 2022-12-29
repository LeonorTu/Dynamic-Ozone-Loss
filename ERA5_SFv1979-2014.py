# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 21:38:56 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/ERA5/o3-1979-2014-JJA-ERA5_2.5.nc')
Va=xr.open_dataset('D:/bylw/ERA5/vwnd-1979-2014-JJA-ERA5_2.5.nc')
#times=O3.time.data[1979:2014]
levs=O3.level.data[:15]  #UTLS区：250hpa-50hpa
lats=O3.latitude.data[:37] #取北半球纬度
a=O3['o3'].loc[:,:,:0,:]  #选取1979-2014共35年的O3数据
# a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
b=Va['v'].loc[:,:,:0,:]  #选取1979-2014共35年的Ua数据
# b1=b.loc[b.time.dt.season=='JJA']  #选取夏季数据

""" 计算STv1 """
a10=np.nanmean(a,axis=3)    #对O3求纬圈平均
a11=np.nanmean(a10,axis=0)   #对O3纬圈平均求时间平均
b10=np.nanmean(b,axis=3)    #对Va求纬圈平均    
b11=np.nanmean(b10,axis=0)   #对Va纬圈平均求时间平均
stv1=np.empty((37,37))  #求STv1
stv1=a11*b11  
stv10=np.empty((37,37))       #求经向散度STv10
dy=6378388.0*2.5*np.pi/180
for x in range(1,36):
    stv10[:,0]=(stv1[:,1]-stv1[:,0])/dy
    stv10[:,x]=(stv1[:,x+1]-stv1[:,x-1])/(2*dy)   
    stv10[:,36]=(stv1[:,36]-stv1[:,35])/dy
stv11=stv10[:,16:24]
# stv12=xr.DataArray(stv11,dims=('z', 'x'))
# stv13=stv12.rolling(z=3, center=True, min_periods=2).mean()  #求滑动平均
stv12=np.nanmean(stv11,axis=1)  

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
plt.hlines(50,-320,350,linewidth=1,color='k',linestyle=':')      
plt.plot(stv12[:15]*(1e+14),levs,linewidth=1,color='#EE82EE',label='SFv1',marker='.')
n.invert_yaxis()          
n.set_yscale('symlog')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
n.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=9.5)
plt.yticks(yticks,ylabel)  
plt.xticks(np.arange(-100,100,20))         
plt.xlim(-100,90)    
n.set_xticklabels(np.arange(-100,100,20),fontfamily="Times New Roman",fontsize=9.5)  

#%%分割cell     
""" 计算STv2 """
a21=np.empty((108,37,37,144))  #对O3求纬偏分量
for y in range(144):
        a21[:,:,:,y]=a[:,:,:,y]-a10
a22=np.nanmean(a21,axis=0)   #对O3纬偏分量求时间平均   
stv2=np.empty((37,37,144))  #求STv2 
for y in range(144):
    stv2[:,:,y]=a22[:,:,y]*b11[:,:]    
stv20=np.empty((37,37,144))   #求经向散度STv20    
for x in range(1,36):
    stv20[:,0,:]=(stv2[:,1,:]-stv2[:,0,:])/dy
    stv20[:,x,:]=(stv2[:,x+1,:]-stv2[:,x-1,:])/(2*dy)   
    stv20[:,36,:]=(stv2[:,36,:]-stv2[:,35,:])/dy
stv21=stv20[:,16:24,83:107]
stv22=np.nanmean(stv21,axis=(1,2))
# stv23=xr.DataArray(stv22,dims=('z', 'x'))
# stv24=stv23.rolling(z=3, center=True, min_periods=2).mean()  #求滑动平均
# stv25=np.nanmean(stv24,axis=1)
# print(stv25)
plt.plot(stv22[:15]*(1e+14),levs,linewidth=1,color='#FFA500',label='SFv2',marker='.') 

#%%分割cell 
""" 计算STv3 """
b21=np.empty((108,37,37,144))  #对Ua求纬偏分量
for y in range(144):
        b21[:,:,:,y]=b[:,:,:,y]-b10
b22=np.nanmean(b21,axis=0)  #对Va纬偏分量求时间平均
stv3=np.empty((37,37,144))  #求STv3
for y in range(144):
    stv3[:,:,y]=b22[:,:,y]*a11[:,:]   
stv30=np.empty((37,37,144))  #求经向散度STv30     
for x in range(1,36):
    stv30[:,0,:]=(stv3[:,1,:]-stv3[:,0,:])/dy
    stv30[:,x,:]=(stv3[:,x+1,:]-stv3[:,x-1,:])/(2*dy)   
    stv30[:,36,:]=(stv3[:,36,:]-stv3[:,35,:])/dy

stv31=stv30[:,16:24,83:107]
stv32=np.nanmean(stv31,axis=(1,2))
# stv33=xr.DataArray(stv32,dims=('z', 'x'))
# stv34=stv33.rolling(z=3, center=True, min_periods=2).mean()
# stv35=np.nanmean(stv34,axis=1)   
plt.plot(stv32[:15]*(1e+14),levs,linewidth=1,color='#3CB371',label='SFv3',marker='.')        
# print(stv35[13])

#%%分割cell 
""" 计算STv4 """    
stv4=np.empty((37,37,144))  #求STv3       
for y in range(144):
    stv4[:,:,y]=a22[:,:,y]*b22[:,:,y]  
stv40=np.empty((37,37,144))   #求经向散度STv40 
for x in range(1,36):
    stv40[:,0,:]=(stv4[:,1,:]-stv4[:,0,:])/dy
    stv40[:,x,:]=(stv4[:,x+1,:]-stv4[:,x-1,:])/(2*dy)   
    stv40[:,36,:]=(stv4[:,36,:]-stv4[:,35,:])/dy
    
stv41=stv40[:,16:24,83:107]
stv42=np.nanmean(stv41,axis=(1,2))
# stv43=xr.DataArray(stv42,dims=('z', 'x'))
# stv44=stv43.rolling(z=3, center=True, min_periods=2).mean()
# stv45=np.nanmean(stv44,axis=1)      
plt.plot(stv42[:15]*(1e+14),levs,linewidth=1,color='#1E90FF',label='SFv4',marker='.')  

#%%分割cell 
""" 计算STv总 """
stv=stv12+stv22+stv32+stv42
np.savetxt('sfvls76.csv',stv,fmt='%.40f',delimiter=',')
plt.plot(stv[:15]*(1e+14),levs,linewidth=1,color='k',label='SFv-Total',marker='.')
tu=plt.legend(bbox_to_anchor=(0,-0.25,1,0.2), loc='lower center',mode='expand',ncol=5,prop={'family' : 'Times New Roman', 'size': 9.5})
tu.get_frame().set_linewidth(0.0)
plt.rc('font', family='Times New Roman')
n.set(xlabel=r'浓度变化$\mathrm{-D}$${_y}$${_-}$${_S}$${_F}$$\mathrm{/(10}$${^-}$${^1}$${^4}$ $\mathrm{mol\;·\;s}$${^-}$${^1}$$\mathrm{\;·\;mol}$${^-}$${^1}$$\mathrm{)}$',ylabel=r'气压$\mathrm{/hPa}$',title=r'$\mathrm{(a)}$ 定常输送')
plt.show()
print(stv[:15])