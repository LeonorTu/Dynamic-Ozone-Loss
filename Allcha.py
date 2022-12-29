# -*- coding: utf-8 -*-
"""
Created on Mon Jul 19 23:00:20 2021

@author: Leonor Tu
"""

import xarray as xr
import numpy as np 
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import t
stu=np.loadtxt('stuwl.csv',delimiter=',')
Stu=np.loadtxt('Stuls.csv',delimiter=',')
ttu=np.loadtxt('ttuwl.csv',delimiter=',')
Ttu=np.loadtxt('Ttuls.csv',delimiter=',')
stv=np.loadtxt('stvwl.csv',delimiter=',')
Stv=np.loadtxt('Stvls.csv',delimiter=',')
ttv=np.loadtxt('ttvwl.csv',delimiter=',')
Ttv=np.loadtxt('Ttvls.csv',delimiter=',')
O3=xr.open_dataset('D:/bylw/o3/o3_Amon_CESM2-WACCM_ssp370_r1i1p1f1_gn_201501-206412.nc')
Ua=xr.open_dataset('D:/bylw/ua/ua_Amon_CESM2-WACCM_ssp370_r1i1p1f1_gn_201501-206412.nc')
#times=O3.time.data[1979:2014]
levs=O3.plev.data[8:14]   #UTLS区：250hpa-50hpa
lats=O3.lat.data[96:] #取北半球纬度
st=stu+stv
St=Stu+Stv
tstst=[]
for z in range(19):
    print(z)
    Avet=np.mean(St[z,:])
    avet=np.mean(st[z,:])
    Vart=np.var(St[z,:])
    vart=np.var(st[z,:])
    tst=(avet-Avet)/np.sqrt((vart+Vart)/15)
    print(tst)
    print('jisuan',stats.ttest_ind(st[z,:],St[z,:],equal_var=False))
    tstst.append(tst)
st1=np.nanmean(st,axis=1)
St1=np.nanmean(St,axis=1)
stcha=st1[8:14]-St1[8:14]
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
plt.plot([0]*6,levs/100.0,linewidth=1,color='k',linestyle='--') 
# plt.hlines(50,-40,16,linewidth=1,color='k',linestyle=':') 
plt.plot(-stcha*(1e+14),levs/100.0,linewidth=1,color='#EE82EE',label='SF-Total',marker='.')
t_alpha= t.ppf(1-0.05/2,28)
for i in range(8,14):
    if abs(tstst[i])>=t_alpha:
        plt.scatter(-stcha[i-8]*(1e+14),levs[i-8]/100,s=200,c='#EE82EE',alpha=0.3,marker='o') 
n.invert_yaxis()          
n.set_yscale('symlog')
yticks=[50,70,100,150,200,250]
ylabel=[50,70,100,150,200,250]
n.set_yticklabels([50,70,100,150,200,250],fontfamily="Times New Roman",fontsize=9.5)
plt.yticks(yticks,ylabel)
plt.xticks(np.arange(-4,4.1,1))         
plt.xlim(-4,4) 
n.set_xticklabels(np.arange(-4,4.1,1),fontfamily="Times New Roman",fontsize=9.5)   

tt=ttu+ttv
Tt=Ttu+Ttv
tsttt=[]
for z in range(19):
    print(z)
    Avet=np.mean(Tt[z,:])
    avet=np.mean(tt[z,:])
    Vart=np.var(Tt[z,:])
    vart=np.var(tt[z,:])
    ttt=(avet-Avet)/np.sqrt((vart+Vart)/15)
    print(ttt)
    print('jisuan',stats.ttest_ind(tt[z,:],Tt[z,:],equal_var=False))
    tsttt.append(ttt)
tt1=np.nanmean(tt,axis=1)
Tt1=np.nanmean(Tt,axis=1)
ttcha=tt1[8:14]-Tt1[8:14]

plt.plot(-ttcha*(1e+14),levs/100.0,linewidth=1,color='#FFA500',label='TF-Total',marker='.')
for i in range(8,14):
    if abs(tsttt[i])>=t_alpha:
        plt.scatter(-ttcha[i-8]*(1e+14),levs[i-8]/100,s=200,c='#FFA500',alpha=0.3,marker='o') 

#%%分割cell
ut=stu+ttu
Ut=Stu+Ttu
tstu=[]
for z in range(19):
    print(z)
    Avet=np.mean(Ut[z,:])
    avet=np.mean(ut[z,:])
    Vart=np.var(Ut[z,:])
    vart=np.var(ut[z,:])
    tut=(avet-Avet)/np.sqrt((vart+Vart)/15)
    print(tut)
    print('jisuan',stats.ttest_ind(ut[z,:],Ut[z,:],equal_var=False))
    tstu.append(tut)
ut1=np.nanmean(ut,axis=1)
Ut1=np.nanmean(Ut,axis=1)
utcha=ut1[8:14]-Ut1[8:14]
plt.plot(-utcha*(1e+14),levs/100.0,linewidth=1,color='#3CB371',label='u-Total',marker='.')
for i in range(8,14):
    if abs(tstu[i])>=t_alpha:
        plt.scatter(-utcha[i-8]*(1e+14),levs[i-8]/100,s=200,c='#3CB371',alpha=0.3,marker='o') 

#%%分割cell
vt=stv+ttv
Vt=Stv+Ttv
tstv=[]
for z in range(19):
    print(z)
    Avet=np.mean(Vt[z,:])
    avet=np.mean(vt[z,:])
    Vart=np.var(Vt[z,:])
    vart=np.var(vt[z,:])
    tvt=(avet-Avet)/np.sqrt((vart+Vart)/15)
    print(tvt)
    print('jisuan',stats.ttest_ind(vt[z,:],Vt[z,:],equal_var=False))
    tstv.append(tvt)
vt1=np.nanmean(vt,axis=1)
Vt1=np.nanmean(Vt,axis=1)
vtcha=vt1[8:14]-Vt1[8:14]
plt.plot(-vtcha*(1e+14),levs/100.0,linewidth=1,color='#1E90FF',label='v-Total',marker='.')
for i in range(8,14):
    if abs(tstv[i])>=t_alpha:
        plt.scatter(-vtcha[i-8]*(1e+14),levs[i-8]/100,s=200,c='#1E90FF',alpha=0.3,marker='o')
        
#%%分割cell
tot=stu+ttu+stv+ttv
Tot=Stu+Ttu+Stv+Ttv
tsto=[]
for z in range(19):
    print(z)
    Avet=np.mean(Tot[z,:])
    avet=np.mean(tot[z,:])
    Vart=np.var(Tot[z,:])
    vart=np.var(tot[z,:])
    to=(avet-Avet)/np.sqrt((vart+Vart)/15)
    print(to)
    print('jisuan',stats.ttest_ind(tot[z,:],Tot[z,:],equal_var=False))
    tsto.append(to)
tot1=np.nanmean(tot,axis=1)
Tot1=np.nanmean(Tot,axis=1)
totcha=tot1[8:14]-Tot1[8:14]
plt.plot(-totcha*(1e+14),levs/100.0,linewidth=1,color='k',label='All',marker='.')
for i in range(8,14):
    if abs(tsto[i])>=t_alpha:
        plt.scatter(-totcha[i-8]*(1e+14),levs[i-8]/100,s=200,c='k',alpha=0.3,marker='o')       
#%%分割cell
tu=plt.legend(bbox_to_anchor=(0,-0.25,1,0.2),loc='lower center',mode='expand',ncol=5,prop={'family' : 'Times New Roman', 'size': 9.5})
tu.get_frame().set_linewidth(0.0)
plt.rc('font', family='Times New Roman')
n.set(xlabel=r'浓度变化差值$\mathrm{-△D}$$\mathrm{/(10}$${^-}$${^1}$${^4}$ $\mathrm{mol\;·\;s}$${^-}$${^1}$$\mathrm{\;·\;mol}$${^-}$${^1}$$\mathrm{)}$',ylabel=r'气压$\mathrm{/hPa}$')
plt.show()