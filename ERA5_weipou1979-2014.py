# -*- coding: utf-8 -*-
"""
Created on Thu Jul  8 10:41:00 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/ERA5/o3-1979-2014-JJA-ERA5_2.5.nc')
levs=O3.level.data[:15]  #UTLS区：250hpa-50hpa
lats=O3.latitude.data[:37] #取北半球纬度
a=O3['o3'].loc[:,:200,:0,:]  #选取1979-2014共35年 250hpa-50hpa的O3数据
# a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
a10=np.nanmean(a,axis=3)  #对O3求纬圈平均
a11=np.empty((108,15,37,144))  #对O3求纬偏分量
for y in range(144):
    a11[:,:,:,y]=a[:,:,:,y]-a10
a12=a11[:,:,:,83:107]
a13=np.nanmean(a12,axis=(0,3))
config = {
    "font.family":'serif',
    "font.size": 9,
    "mathtext.fontset":'stix',
    "font.serif": ['STZhongsong'],
}
plt.rcParams.update(config)
plt.rcParams['axes.unicode_minus']=False
plt.rcParams['figure.figsize'] = (3.3,2.23) # 单位是inches
fig=plt.figure(dpi=800)  #添加画布
ax=fig.add_axes([0,0,1,1])  #添加子图
plt.hlines(50,0,180,linewidth=1,color='k',linestyle=':')
ax.invert_yaxis()  #反转纵轴
ax.set_yscale('symlog')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
plt.yticks(yticks,ylabel) 
ax.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=8)
plt.xticks(np.arange(0,91,10))
plt.xlim(0,90)
ax.set_xticklabels([r'0$^\degree$',r'10$^\degree$N',r'20$^\degree$N',r'30$^\degree$N',r'40$^\degree$N',
                  r'50$^\degree$N',r'60$^\degree$N',r'70$^\degree$N',r'80$^\degree$N',r'90$^\degree$N'],fontfamily="Times New Roman",fontsize=8)#转换为纬度格式
XX,YY=np.meshgrid(lats,levs,indexing='ij')
ca=ax.contourf(XX,YY,a13.T*(1e+8),levels=np.arange(-45,45.1,5),cmap=plt.cm.RdBu_r)
plt.rc('font', family='Times New Roman')
cb=plt.colorbar(ca,orientation='horizontal',fraction=0.05,pad=0.1)
cb.ax.tick_params(labelsize=8)
ax.set(ylabel=r'气压$\mathrm{/hPa}$',title='$\mathrm{(b)}$ 纬向垂直分布')
plt.show()
