# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 11:54:56 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/cmip6_3/o3_Amon_CMIP6_3models_historical_r1i1p1f1_1979-2014-JJA-2.5.nc')
Ua=xr.open_dataset('D:/bylw/cmip6_3/ua_Amon_CMIP6_3models_historical_r1i1p1f1_1979-2014-JJA-2.5.nc')
#times=O3.time.data[1979:2014]
levs=O3.plev.data[9:]   #UTLS区：250hpa-50hpa
lats=O3.lat.data[36:] #取北半球纬度
a=O3['o3'].loc[:,:,0:,:]  #选取1979-2014共35年北半球的O3数据
# a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
b=Ua['ua'].loc[:,:,0:,:]  #选取1979-2014共35年北半球的Ua数据
# b1=b.loc[b.time.dt.season=='JJA']  #选取夏季数据
#%%分割cell
""" 计算STu1 """
a10=np.nanmean(a,axis=3)    #对O3求纬圈平均
a11=np.nanmean(a10,axis=0)   #对O3纬圈平均求时间平均
b10=np.nanmean(b,axis=3)    #对Ua求纬圈平均    
b11=np.nanmean(b10,axis=0)   #对Ua纬圈平均求时间平均
stu1=np.empty((19,37))  #求STu1
stu1=a11*b11  

config = {
    "font.family":'serif',
    "font.size": 10.5,
    "mathtext.fontset":'stix',
    "font.serif": ['STZhongsong'],
}
plt.rcParams.update(config)
plt.rcParams['axes.unicode_minus']=False
fig=plt.figure(1,dpi=800)   #画图
plt.rcParams['figure.figsize'] = (4,3) # 单位是inches
n=fig.add_axes([0, 0, 1, 1])     
plt.plot([0]*10,levs/100.0,linewidth=1,color='k',linestyle='--') 
plt.hlines(50,-150,20,linewidth=1,color='k',linestyle=':')      
plt.plot([0]*10,levs/100.0,linewidth=1,color='#EE82EE',label='SFu1',marker='.')
n.invert_yaxis()          
n.set_yscale('symlog')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
n.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=9.5)
plt.yticks(yticks,ylabel) 
plt.xticks(np.arange(-130,21,10))     
plt.xlim(-130,20)
n.set_xticklabels(np.arange(-130,21,10),fontfamily="Times New Roman",fontsize=9.5)       

#%%分割cell 
""" 计算STu2 """
a21=np.empty((108,19,37,144))  #对O3求纬偏分量
for y in range(144):
        a21[:,:,:,y]=a[:,:,:,y]-a10
a22=np.nanmean(a21,axis=0)   #对O3纬偏分量求时间平均   
stu2=np.empty((19,37,144))  #求STu2 
for y in range(144):
    stu2[:,:,y]=a22[:,:,y]*b11[:,:]    
stu20=np.empty((19,37,144))   #求纬向散度STu20    
dlon=2.5*np.pi/180.0
for x in range(37):  
    dx=6378388.0*np.cos(lats[x]*(np.pi/180.0))*dlon
    for y in range(1,143):
        stu20[:,x,0]=(stu2[:,x,1]-stu2[:,x,0])/dx
        stu20[:,x,y]=(stu2[:,x,y+1]-stu2[:,x,y-1])/(2*dx)  
        stu20[:,x,143]=(stu2[:,x,143]-stu2[:,x,142])/dx

stu21=stu20[:,12:19,12:40]
stu22=np.nanmean(stu21,axis=(1,2))
# stu23=xr.DataArray(stu22,dims=('z', 'x'))
# stu24=stu23.rolling(z=3, center=True, min_periods=2).mean()  #对高度求三点滑动平均
# stu25=np.nanmean(stu24,axis=1) 
plt.plot(-stu22[9:]*(1e+14),levs/100.0,linewidth=1,color='#FFA500',label='SFu2',marker='.') 

#%%分割cell 
""" 计算STu3 """
b21=np.empty((108,19,37,144))  #对Ua求纬偏分量
for y in range(144):
    b21[:,:,:,y]=b[:,:,:,y]-b10
b22=np.nanmean(b21,axis=0)  #对Ua纬偏分量求时间平均
stu3=np.empty((19,37,144))  #求STu3
for y in range(144):
    stu3[:,:,y]=b22[:,:,y]*a11[:,:]   
stu30=np.empty((19,37,144))  #求纬向散度STu30     
for x in range(37):  
    dx=6378388.0*np.cos(lats[x]*(np.pi/180.0))*dlon
    for y in range(1,143):
        stu30[:,x,0]=(stu3[:,x,1]-stu3[:,x,0])/dx
        stu30[:,x,y]=(stu3[:,x,y+1]-stu3[:,x,y-1])/(2*dx)  
        stu30[:,x,143]=(stu3[:,x,143]-stu3[:,x,142])/dx

stu31=stu30[:,12:19,12:40]
stu32=np.nanmean(stu31,axis=(1,2))  
# stu33=xr.DataArray(stu32,dims=('z', 'x'))
# stu34=stu33.rolling(z=3, center=True, min_periods=2).mean()  #对高度求三点滑动平均
# stu35=np.nanmean(stu34,axis=1)      
plt.plot(-stu32[9:]*(1e+14),levs/100.0,linewidth=1,color='#3CB371',label='SFu3',marker='.')        
# print(stu35[15]) 

#%%分割cell 
""" 计算STu4 """    
stu4=np.empty((19,37,144))  #求STu3       
for y in range(144):
    stu4[:,:,y]=a22[:,:,y]*b22[:,:,y]  
stu40=np.empty((19,37,144))   #求纬向散度STu40 
for x in range(37):  
    dx=6378388.0*np.cos(lats[x]*(np.pi/180.0))*dlon
    for y in range(1,143):
        stu40[:,x,0]=(stu4[:,x,1]-stu4[:,x,0])/dx
        stu40[:,x,y]=(stu4[:,x,y+1]-stu4[:,x,y-1])/(2*dx)  
        stu40[:,x,143]=(stu4[:,x,143]-stu4[:,x,142])/dx   

stu41=stu40[:,12:19,12:40]
stu42=np.nanmean(stu41,axis=(1,2))   
# stu43=xr.DataArray(stu42,dims=('z', 'x'))
# stu44=stu43.rolling(z=3, center=True, min_periods=2).mean()  #对高度求三点滑动平均
# stu45=np.nanmean(stu44,axis=1)             
plt.plot(-stu42[9:]*(1e+14),levs/100.0,linewidth=1,color='#1E90FF',label='SFu4',marker='.')  

""" 计算STu总 """
stu=stu22+stu32+stu42
np.savetxt('sfuls74.csv',stu,fmt='%.40f',delimiter=',')
plt.plot(-stu[9:]*(1e+14),levs/100.0,linewidth=1,color='k',label='SFu-Total',marker='.')
tu=plt.legend(bbox_to_anchor=(0,-0.25,1,0.2), loc='lower center',mode='expand',ncol=5,prop={'family' : 'Times New Roman', 'size': 9.5})
tu.get_frame().set_linewidth(0.0)
plt.rc('font', family='Times New Roman')
n.set(xlabel=r'浓度变化$\mathrm{-D}$${_x}$${_-}$${_S}$${_F}$$\mathrm{/(10}$${^-}$${^1}$${^4}$ $\mathrm{mol\;·\;s}$${^-}$${^1}$$\mathrm{\;·\;mol}$${^-}$${^1}$$\mathrm{)}$',ylabel=r'气压$\mathrm{/hPa}$',title=r'$\mathrm{(a)}$ 定常输送')
plt.show()
print(stu[9:])