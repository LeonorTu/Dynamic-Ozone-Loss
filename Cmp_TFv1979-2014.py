# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 21:58:26 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/cmip6_3/o3_Amon_CMIP6_3models_historical_r1i1p1f1_1979-2014-JJA-2.5.nc')
Va=xr.open_dataset('D:/bylw/cmip6_3/va_Amon_CMIP6_3models_historical_r1i1p1f1_1979-2014-JJA-2.5.nc')
#times=O3.time.data[1979:2014]
levs=O3.plev.data[9:]   #UTLS区：250hpa-50hpa
lats=O3.lat.data[36:] #取北半球纬度
a=O3['o3'].loc[:,:,0.:,:]  #选取1979-2014共35年的O3数据
# a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
b=Va['va'].loc[:,:,0.:,:]  #选取1979-2014共35年的Ua数据
# b1=b.loc[b.time.dt.season=='JJA']  #选取夏季数据

""" 计算TTv1 """
a10=np.nanmean(a,axis=0)    #对O3求时间平均
a11=np.empty((108,19,37,144))  #对O3求距平分量
for t in range(108):
    a11[t,:,:,:]=a[t,:,:,:]-a10
a20=np.nanmean(a11,axis=3)  #对O3距平分量求纬圈平均
b10=np.nanmean(b,axis=0)    #对Va求时间平均
b11=np.empty((108,19,37,144))  #对Va求距平分量
for t in range(108):
    b11[t,:,:,:]=b[t,:,:,:]-b10
b20=np.nanmean(b11,axis=3)  #对Va距平分量求纬圈平均
ttv1=np.empty((108,19,37))  #求TTv1
ttv1=a20*b20
ttv10=np.nanmean(ttv1,axis=0)  #对ttv1求时间平均
ttv11=np.empty((19,37))  #求经向散度TTv11
dy=6378388.0*2.5*np.pi/180
for x in range(1,36):
    ttv11[:,0]=(ttv10[:,1]-ttv10[:,0])/dy
    ttv11[:,x]=(ttv10[:,x+1]-ttv10[:,x-1])/(2*dy)
    ttv11[:,36]=(ttv10[:,36]-ttv10[:,35])/dy   
ttv12=ttv11[:,12:19]
# ttv13=xr.DataArray(ttv12,dims=('z', 'x'))
# ttv14=ttv13.rolling(z=3, center=True, min_periods=2).mean()  #求滑动平均
ttv13=np.nanmean(ttv12,axis=1)
# print(ttv15)

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
plt.plot([0]*10,levs/100.0,linewidth=1,color='k',linestyle='--')    
plt.hlines(50,-40,16,linewidth=1,color='k',linestyle=':')     
plt.plot(-ttv13[9:]*(1e+14),levs/100.0,linewidth=1,color='#EE82EE',label='TFv1',marker='.') 
n.invert_yaxis()          
n.set_yscale('symlog')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
n.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=9.5)
plt.yticks(yticks,ylabel)
plt.xticks(np.arange(-1.2,1.81,0.2))         
plt.xlim(-1.2,1.81)  
n.set_xticklabels([-1.2,-1.0,-0.8,-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8],fontfamily="Times New Roman",fontsize=9.5)  
 
#%%分割cell        
""" 计算TTv2 """      
b21=np.empty((108,19,37,144))  #对Va距平分量求纬圈距平
for y in range(144):
    b21[:,:,:,y]=b11[:,:,:,y]-b20
ttv2=np.empty((108,19,37,144))  #求TTv2 
for y in range(144):
    ttv2[:,:,:,y]=b21[:,:,:,y]*a20  
ttv20=np.nanmean(ttv2,axis=0)    #对TTv2求时间平均
ttv21=np.empty((19,37,144))   #求纬向散度TTv20    
for x in range(1,36):
    ttv21[:,0,:]=(ttv20[:,1,:]-ttv20[:,0,:])/dy
    ttv21[:,x,:]=(ttv20[:,x+1,:]-ttv20[:,x-1,:])/(2*dy) 
    ttv21[:,36,:]=(ttv20[:,36,:]-ttv20[:,35,:])/dy

ttv22=ttv21[:,12:19,16:37]
ttv23=np.nanmean(ttv22,axis=(1,2))
# ttv24=xr.DataArray(ttv23,dims=('z', 'x'))
# ttv25=ttv24.rolling(z=3, center=True, min_periods=2).mean()  #求滑动平均
# ttv26=np.nanmean(ttv25,axis=1)       
plt.plot(-ttv23[9:]*(1e+14),levs/100.0,linewidth=1,color='#FFA500',label='TFv2',marker='.') 

#%%分割cell  
""" 计算TTv3 """
a21=np.empty((108,19,37,144))  #对O3距平分量求纬圈距平
for y in range(144):
    a21[:,:,:,y]=a11[:,:,:,y]-a20
ttv3=np.empty((108,19,37,144))  #求TTv3 
for y in range(144):
    ttv3[:,:,:,y]=a21[:,:,:,y]*b20  
ttv30=np.nanmean(ttv3,axis=0)    #对TTv2求时间平均
ttv31=np.empty((19,37,144))   #求纬向散度TTv31      
for x in range(1,36):
    ttv31[:,0,:]=(ttv30[:,1,:]-ttv30[:,0,:])/dy
    ttv31[:,x,:]=(ttv30[:,x+1,:]-ttv30[:,x-1,:])/(2*dy) 
    ttv31[:,36,:]=(ttv30[:,36,:]-ttv30[:,35,:])/dy

ttv32=ttv31[:,12:19,16:37]
ttv33=np.nanmean(ttv32,axis=(1,2))
# ttv34=xr.DataArray(ttv33,dims=('z', 'x'))
# ttv35=ttv34.rolling(z=3, center=True, min_periods=2).mean()
# ttv36=np.nanmean(ttv35,axis=1)
plt.plot(-ttv33[9:]*(1e+14),levs/100.0,linewidth=1,color='#3CB371',label='TFv3',marker='.')         

#%%分割cell  
""" 计算TTv4 """  
ttv4=np.empty((108,19,37,144))  #求TTv4 
for y in range(144):
    ttv4[:,:,:,y]=a21[:,:,:,y]*b21[:,:,:,y]   
ttv40=np.nanmean(ttv4,axis=0)    #对TTv2求时间平均
ttv41=np.empty((19,37,144))   #求纬向散度TTv20        
for x in range(1,36):
    ttv41[:,0,:]=(ttv40[:,1,:]-ttv40[:,0,:])/dy
    ttv41[:,x,:]=(ttv40[:,x+1,:]-ttv40[:,x-1,:])/(2*dy) 
    ttv41[:,36,:]=(ttv40[:,36,:]-ttv40[:,35,:])/dy
ttv42=ttv41[:,12:19,16:37]
ttv43=np.nanmean(ttv42,axis=(1,2))
# ttv44=xr.DataArray(ttv43,dims=('z', 'x'))
# ttv45=ttv44.rolling(z=3, center=True, min_periods=2).mean()
# ttv46=np.nanmean(ttv45,axis=1)
plt.plot(-ttv43[9:]*(1e+14),levs/100.0,linewidth=1,color='#1E90FF',label='TFv4',marker='.') 

#%%分割cell 
""" 计算TTv总 """
ttv=ttv13+ttv23+ttv33+ttv43
np.savetxt('tfvls74.csv',ttv,fmt='%.40f',delimiter=',')
plt.plot(-ttv[9:]*(1e+14),levs/100.0,linewidth=1,color='k',label='TFv-Total',marker='.')
tu=plt.legend(bbox_to_anchor=(0,-0.25,1,0.2), loc='lower center',mode='expand',ncol=5,prop={'family' : 'Times New Roman', 'size': 9.5})
tu.get_frame().set_linewidth(0.0)
plt.rc('font', family='Times New Roman')
n.set(xlabel=r'浓度变化$\mathrm{-D}$${_y}$${_-}$${_T}$${_F}$$\mathrm{/(10}$${^-}$${^1}$${^4}$ $\mathrm{mol\;·\;s}$${^-}$${^1}$$\mathrm{\;·\;mol}$${^-}$${^1}$$\mathrm{)}$',ylabel=r'气压$\mathrm{/hPa}$',title=r'$\mathrm{(b)}$ 瞬变输送')
plt.show()
print(ttv[9:])