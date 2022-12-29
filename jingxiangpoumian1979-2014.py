# -*- coding: utf-8 -*-
"""
Created on Sat May 15 16:02:58 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
O3=xr.open_dataset('D:/bylw/o3/o3_Amon_CESM2-WACCM_historical_r1i1p1f1_gn_185001-201412.nc')
levs=O3.plev.data[9:]   #250hpa-
lats=O3.lat.data[96:192] #取北半球纬度
lons=O3.lon.data[:145]
a=O3['o3'].loc['1979-01-15':'2014-12-15',20000.:,0.:,:]  #选取1979-2014共35年 250hpa-50hpa的O3数据
a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
a10=np.nanmean(a1,axis=3)  #对O3求纬圈平均
a11=np.empty((108,10,96,288))  #对O3求纬偏分量
for y in range(288):
    a11[:,:,:,y]=a1[:,:,:,y]-a10
a12=a11[:,:,26:49,:145]
a13=np.nanmean(a12,axis=(0,2))
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
ax.set_yscale('log')
yticks=[1,5,10,20,30,50,70,100,150,200]
ylabel=[1,5,10,20,30,50,70,100,150,200]
plt.yticks(yticks,ylabel) 
ax.set_yticklabels([1,5,10,20,30,50,70,100,150,200],fontfamily="Times New Roman",fontsize=8)
plt.xticks(np.arange(0,181,20))
plt.xlim(0,180)
ax.set_xticklabels([r'0$^\degree$',r'20$^\degree$E',r'40$^\degree$E',r'60$^\degree$E',r'80$^\degree$E',
                 r'100$^\degree$E',r'120$^\degree$E',r'140$^\degree$E',r'160$^\degree$E',r'180$^\degree$'],fontfamily="Times New Roman",fontsize=8)#转换为经度格式
XX,YY=np.meshgrid(lons,levs/100,indexing='ij')
ca=ax.contourf(XX,YY,a13.T*(1e+8),levels=np.arange(-32,32.1,4),cmap=plt.cm.RdBu_r)
plt.rc('font', family='Times New Roman')
cb=plt.colorbar(ca,orientation='horizontal',fraction=0.05,pad=0.1)
cb.ax.tick_params(labelsize=8)
ax.set(ylabel=r'气压$\mathrm{/hPa}$',title='$\mathrm{(c)}$ 经向垂直分布')
plt.show()