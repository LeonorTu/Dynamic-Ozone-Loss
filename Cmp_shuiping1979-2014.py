# -*- coding: utf-8 -*-
"""
Created on Sat Jul  3 19:29:47 2021

@author: Leonor Tu
"""

import cartopy.crs  as ccrs
import matplotlib.pyplot as plt
import numpy as  np
import xarray as xr  
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter,LatitudeFormatter
import cartopy.io.shapereader as shpreader  
from matplotlib.font_manager import FontProperties
O3=xr.open_dataset('D:/bylw/cmip6_3/o3_Amon_CMIP6_3models_historical_r1i1p1f1_1979-2014-JJA-2.5.nc')
levs=O3.plev.data[9:]   #250hpa-
lats=O3.lat.data[36:] #取北半球纬度
lons=O3.lon.data 
a=O3['o3'].loc[:,20000.:,0.:,:]  #选取1979-2014共35年 200hpa- 北半球的O3数据
# a1=a.loc[a.time.dt.season=='JJA']  #选取夏季数据
a10=np.nanmean(a,axis=3)  #对O3求纬圈平均
a11=np.empty((108,10,37,144))  #对O3求纬偏分量
for y in range(144):
    a11[:,:,:,y]=a[:,:,:,y]-a10
a12=np.nanmean(a11,axis=(0,1))
config = {
    "font.family":'serif',
    "font.size": 9,
    "mathtext.fontset":'stix',
    "font.serif": ['STZhongsong'],
}
plt.rcParams.update(config)
plt.rcParams['axes.unicode_minus']=False
plt.rcParams['figure.figsize'] = (3.3,2.5) # 单位是inches
fig=plt.figure(dpi=800)  #添加画布
leftlon, rightlon, lowerlat, upperlat = (0,180,15,55)
ax = fig.add_axes([0,0,1,1],projection=ccrs.PlateCarree())
ax.set_extent([leftlon, rightlon, lowerlat, upperlat], crs=ccrs.PlateCarree())
ax.add_feature(cfeature.COASTLINE.with_scale('110m'))
ax.add_feature(cfeature.LAKES, alpha=0.5)
ax.set_xticks(np.arange(leftlon,rightlon+10,10), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lowerlat,upperlat+3,5), crs=ccrs.PlateCarree())
lianghe = shpreader.Reader('D:/bylw/TP/DBATP_Line.shp').geometries()
ax.add_geometries(lianghe, ccrs.PlateCarree(),facecolor ='none',edgecolor ='sienna',lw=0.8,zorder=1,alpha=6)
XX,YY=np.meshgrid(lons[:72],lats,indexing='ij')
ca=ax.contourf(XX,YY,(a12[:,:72]).T*(1e+8),levels=np.arange(-14,14.1,2),cmap=plt.cm.RdBu_r)
plt.xticks(np.arange(0,181,20))
plt.xlim(0,180)
ax.set_xticklabels([r'0$^\degree$',r'20$^\degree$E',r'40$^\degree$E',r'60$^\degree$E',r'80$^\degree$E',r'100$^\degree$E',r'120$^\degree$E',r'140$^\degree$E',r'160$^\degree$E',r'180$^\degree$'],fontfamily="Times New Roman",fontsize=8)
plt.yticks(np.arange(0,91,10))
ax.set_yticklabels([r'0$^\degree$',r'10$^\degree$N',r'20$^\degree$N',r'30$^\degree$N',r'40$^\degree$N',r'50$^\degree$N',r'60$^\degree$N',r'70$^\degree$N',r'80$^\degree$N',r'90$^\degree$N'],fontfamily="Times New Roman",fontsize=8)
plt.rc('font', family='Times New Roman')
cb=plt.colorbar(ca,orientation='horizontal',fraction=0.05,pad=0.12)
cb.ax.tick_params(labelsize=8)
ax.set(title=r'$\mathrm{(a)}$ 水平分布')
plt.show()
