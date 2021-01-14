#%%
# Imports
#################################################################
from netCDF4 import Dataset
import numpy as np
import numpy.ma as ma
from scipy.optimize import curve_fit
import statsmodels.api as sm
from statsmodels.tsa.seasonal import seasonal_decompose
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import cartopy.crs as ccrs
#################################################################

#%%
time = []

thetaF = np.zeros((173, 360, 42, 828))
T_error = np.zeros((173, 360, 42, 828))

k=0

#data_dir = r'C:\Users\jack3\Documents\Imperial College London\Third Year\BSc Project - Ocean Heat Uptake\Data\\'
data_dir = r'/rds/general/user/jc9017/home/work//' #need the r at the front otherwise struggles with \
for n in range(1950, 2019):
    for m in range(1,13):
        if m < 10:
            filename = 'EN.4.2.1.f.analysis.g10.{}0{}.nc'.format(n,m)
            print(filename)
        else:
            filename = 'EN.4.2.1.f.analysis.g10.{}{}.nc'.format(n,m)
            print(filename)

        dataset = Dataset(data_dir + filename)
        x = dataset.variables['lon'][:]
        y = dataset.variables['lat'][:]
        z = dataset.variables['depth'][:]
        t = dataset.variables['time'][:]
        theta = dataset.variables['temperature'][:] # shape: theta(t,z,y,x)
        sigma_t = dataset.variables['temperature_uncertainty'][:] # shape: sigma_t(t,z,y,x)
        
        
        time.append(t[0])
        
	for i in range(360):
        	for j in range(173):
                	for p in range(42):
                    thetaF[j,i,p,k] = theta[0,p,j,i]   
                    T_error[j,i,p,k] = sigma_t[0,p,j,i]
        
        k += 1

#%%

# Plot 1: Overall Time Series

HC = np.zeros((173,360,828))
HC_global = []
HC_middle = []
HC_upper = []
HC_lower = []
factor = (np.pi/180)**2

for t in range(828):
    heat = []
    heat_u = []
    heat_l = []
    heat_m = []
    for i in range(360):
        for j in range(173):
            j_ = y[j]*((np.pi)/180)
            A = (factor*(6730e3)**2)*np.cos(j_)*np.trapz(thetaF[j,i,:,t],z[:])
            m = (factor*(6730e3)**2)*np.cos(j_)*np.trapz(thetaF[j,i,:25,t],z[:25])
            u = (factor*(6730e3)**2)*np.cos(j_)*np.trapz((thetaF[j,i,:25,t]+T_error[j,i,:25,t]),z[:25])
            l = (factor*(6730e3)**2)*np.cos(j_)*np.trapz((thetaF[j,i,:25,t]-T_error[j,i,:25,t]),z[:25])       
            HC[j,i,t] = 1025*3850*A
            heat.append(1025*3850*A)
            heat_u.append(1025*3850*u)
            heat_l.append(1025*3850*l)
            heat_m.append(1025*3850*m)
    
    heat_ = ma.masked_invalid(heat)
    heat_u_ = ma.masked_invalid(heat_u)
    heat_l_ = ma.masked_invalid(heat_l)
    heat_m_ = ma.masked_invalid(heat_m)
    HC_global.append(heat_.sum())
    HC_upper.append(heat_u_.sum())
    HC_lower.append(heat_l_.sum())
    HC_middle.append(heat_m_.sum())

#%%
#2 cells producing value for total heat flux with 2 different methods    
def linear_fit(x,a,b):
    y = a*x + b
    return y
time__ = np.linspace(1,828,828)
HC_p0 = (1e20,1e28)
HCt_param, HCt_cov = curve_fit(linear_fit, time, HC_global,p0=HC_p0)
HCt_regress = HCt_param[0]/((86400)*4*np.pi*(6370e3)**2)
print('Total Flux in Wm^2',HCt_regress)       
print('Total FLux Error',np.sqrt(HCt_cov[0]))
#%%
HCt_anomaly = HC_global - HC_global[0]
HCt_anom_middle = HC_middle - np.mean(HC_middle)
HCt_anom_upper = HC_upper - np.mean(HC_upper)
HCt_anom_lower = HC_lower - np.mean(HC_lower)

#%%
time_ = np.arange('1950', '2019', dtype='datetime64[M]')
  
fig, ax = plt.subplots()
fig.set_size_inches(10,5)
ax.plot(time_, (HCt_anomaly), color = 'mediumblue', lw=3, alpha=0.6)
ax.axvline(x='2000-01', color='k', alpha=0.4, lw=2, ls='--')
#ax.axvline(x='1902-10', color='k', alpha=0.4, lw=2, ls='--')
#ax.axvline(x='1912-06', color='k', alpha=0.4, lw=2, ls='--')
ax.axvline(x='1991-06', color='k', alpha=0.4, lw=2, ls='--')
ax.axvline(x='1982-04', color='k', alpha=0.4, lw=2, ls='--')
ax.axvline(x='1963-02', color='k', alpha=0.4, lw=2, ls='--')
ax.axvline(x='1997-07', color='k', alpha=0.4, lw=2, ls='--')
ax.axvline(x='1982-02', color='k', alpha=0.4, lw=2, ls='--')
#add other v-lines
#ax.text('1902-12',1.5e23,'Santa Maria; Guatemala',rotation=90)
#ax.text('1912-08',1.5e23,'Novarupta; Alaska',rotation=90)
ax.text('1991-08',1.5e23,'Pinatubo; Philippines',rotation=90)
ax.text('1982-06',1.5e23,'El Chichon; Mexico',rotation=90)
ax.text('1963-04',1.5e23,'Mt Agung; Indonesia',rotation=90)
ax.text('1997-09',1.5e23,'Large El Niño Event',rotation=90)
ax.text('1982-04',1.5e23,'Large El Niño Event',rotation=90)
ax.text('2000-03',1.5e23,'Start of ARGO Data',rotation=90)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)
plt.minorticks_on()
ax.set_ylabel('Ocean Heat Content (Joules)', fontsize=15)
ax.set_xlabel('Time in Years', fontsize=15)
plt.savefig('/rds/general/user/jc9017/home/Total_Heat_Content_PlotNEW.png')
#%%
## Plot 1.5: 1000m Error
fig, ax = plt.subplots()
fig.set_size_inches(10,5)
ax.plot(time_, HCt_anom_upper, color='cornflowerblue', alpha=0.15)
ax.plot(time_, HCt_anom_lower, color='cornflowerblue', alpha=0.15)
ax.fill_between(time_, HCt_anom_upper, HCt_anom_lower, color='cornflowerblue', alpha=0.15)
decomp = sm.tsa.seasonal_decompose(HCt_anom_middle, model = 'additive', freq=12)
ax.plot(time_, decomp, color = 'mediumblue', lw=3, alpha=0.6)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.yticks(fontsize = 12)
plt.xticks(fontsize = 12)
plt.minorticks_on()
ax.set_ylabel('Ocean Heat Content Anomaly (Joules)', fontsize=15)
ax.set_xlabel('Time in Years', fontsize=15)
plt.savefig('/rds/general/user/jc9017/home/1000m_Error_PlotNEW.png')
#%%
## Plot 2: HC % Change v Time for each Layer

HC_global_L1 = []
HC_global_L2 = []
HC_global_L3 = []
HC1 = np.zeros((173,360,828))
HC2 = np.zeros((173,360,828))
HC3 = np.zeros((173,360,828))
tick = 0
for t in range(828):
    tick += 1
    heat1 = []
    heat2 = []
    heat3 = []
    for i in range(360):
        for j in range(173):
            j_ = y[j]*((np.pi)/180)
            A1 = ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)*np.trapz(thetaF[j,i,0:24,t],z[0:24])
            A2 = ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)*np.trapz(thetaF[j,i,24:30,t],z[24:30])
            A3 = ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)*np.trapz(thetaF[j,i,30:42,t],z[30:42])
            HC1[j,i,t] = 1025*3850*A1
            HC2[j,i,t] = 1025*3850*A2
            HC3[j,i,t] = 1025*3850*A3
            heat1.append(1025*3850*A1)
            heat2.append(1025*3850*A2)
            heat3.append(1025*3850*A3)
    
    heat1_ = ma.masked_invalid(heat1)
    heat2_ = ma.masked_invalid(heat2)
    heat3_ = ma.masked_invalid(heat3)
    HC_global_L1.append(heat1_.sum())
    HC_global_L2.append(heat2_.sum())
    HC_global_L3.append(heat3_.sum())
    

#%%
HC_change_L1 = np.zeros(828)
HC_change_L2 = np.zeros(828)
HC_change_L3 = np.zeros(828)
#%%
HC_trend_L1 = sm.tsa.seasonal_decompose(HC_global_L1, model='additive', freq=12)
HC_trend_L1 = ma.masked_invalid(HC_trend_L1.trend)
HC_trend_L2 = sm.tsa.seasonal_decompose(HC_global_L2, model='additive', freq=12)
HC_trend_L2 = ma.masked_invalid(HC_trend_L2.trend)
HC_trend_L3 = sm.tsa.seasonal_decompose(HC_global_L3, model='additive', freq=12)
HC_trend_L3 = ma.masked_invalid(HC_trend_L3.trend)

#%%

for i in range(828):
    HC_change_L1[i] = (HC_trend_L1[i] - ma.mean(HC_trend_L1))
    HC_change_L2[i] = (HC_trend_L2[i] - ma.mean(HC_trend_L2))
    HC_change_L3[i] = (HC_trend_L3[i] - ma.mean(HC_trend_L3))


#%% 
plt.clf()   
fig, ax = plt.subplots()
fig.set_size_inches(10,5)
ax.plot(time_, HC_change_L1, color='lightskyblue', label='0 - 700m')
ax.plot(time_, HC_change_L2, color='cornflowerblue', label='700 - 2000m')
ax.plot(time_, HC_change_L3, color='royalblue', label='2000m+')
ax.fill_between(time_, HC_change_L1, color='lightskyblue')
ax.fill_between(time_, HC_change_L2, color='cornflowerblue')
ax.fill_between(time_, HC_change_L3, color='royalblue')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.set_xlabel('Time in Years', fontsize=15)
ax.set_ylabel('Total OHC Anomaly (Joules)', fontsize=15)
ax.legend()
plt.savefig('/rds/general/user/jc9017/home/Layer_Heat_Content_PlotNEW.png')

#%%
#Heat content for heat maps

HHC1 = np.zeros((173,360,828))
HHC2 = np.zeros((173,360,828))
HHC3 = np.zeros((173,360,828))
tick = 0
for t in range(828):
    tick += 1
    for i in range(360):
        for j in range(173):
            j_ = y[j]*((np.pi)/180)
            A1 = np.trapz(thetaF[j,i,0:24,t], z[0:24])
            A2 = np.trapz(thetaF[j,i,24:30,t], z[24:30])
            A3 = np.trapz(thetaF[j,i,30:42,t], z[30:42])
            HHC1[j,i,t] = 1025*3850*A1
            HHC2[j,i,t] = 1025*3850*A2
            HHC3[j,i,t] = 1025*3850*A3

#%%
#Plot 5: 3x2 Subplot for Heat Maps 

def linear_fit(x,a,b):
    y = a*x + b
    return y

HC1_ = np.where(np.isnan(HC1), None, HHC1)
HC1__ = ma.masked_where(HC1_ == None, HC1_)
HC2_ = np.where(np.isnan(HC2), None, HHC2)
HC2__ = ma.masked_where(HC2_ == None, HC2_)
HC3_ = np.where(np.isnan(HC3), None, HHC3)
HC3__ = ma.masked_where(HC3_ == None, HC3_)

HC1_regress = np.zeros((173,360))
HC2_regress = np.zeros((173,360))
HC3_regress = np.zeros((173,360))
cov_L1 = np.zeros((173,360))
cov_L2 = np.zeros((173,360))
cov_L3 = np.zeros((173,360))

Aet = (86400)   #Use this in heatmap plot


HC_p0 = (5e21, 1e23)

for i in range(360):
    for j in range(173):
        j_ = y[j]*((np.pi)/180)
        if np.isnan(HC3__[j,i,0]) == False: 
    
            HC1___ = np.divide(HC1__[j,i],Aet)
            HC1_param, HC1_cov = curve_fit(linear_fit, time, HC1___, p0 = HC_p0)
            HC1_regress[j,i] = HC1_param[0]
            cov_L1[j,i] = (np.sqrt(HC1_cov[0,0]))#/HC1_regress[j,i])*100
            
            HC2___ = np.divide(HC2__[j,i],Aet)
            HC2_param, HC2_cov = curve_fit(linear_fit, time, HC2___, p0 = HC_p0)
            HC2_regress[j,i] = HC2_param[0]
            cov_L2[j,i] = (np.sqrt(HC2_cov[0,0]))#/HC2_regress[j,i])*100
            
            HC3___ = np.divide(HC3__[j,i],Aet)
            HC3_param, HC3_cov = curve_fit(linear_fit, time, HC3___, p0 = HC_p0)
            HC3_regress[j,i] = HC3_param[0]
            cov_L3[j,i] = (np.sqrt(HC3_cov[0,0]))#/HC3_regress[j,i])*100
        else:
            HC1_regress[j,i] = None
            cov_L1[j,i] = None
            HC2_regress[j,i] = None
            cov_L2[j,i] = None
            HC3_regress[j,i] = None
            cov_L3[j,i] = None

#%%
plt.clf()
fig, ax = plt.subplots(ncols=2, nrows=3, subplot_kw={'projection': ccrs.PlateCarree()})
fig.set_size_inches(10,5)
cmap1 = plt.cm.seismic
cmap2 = plt.cm.magma

cm00 = ax[0,0].contourf(x, y, HC1_regress, 100, cmap=cmap1, transform=ccrs.PlateCarree())
ax[0,0].coastlines()
cm00.set_clim(-1*np.max(HC1_regress),np.max(HC1_regress))
fig.colorbar(cm00, label='$W/m^2$', ax=ax[0,0])


cm01 = ax[0,1].contourf(x, y, cov_L1, 10, cmap=cmap2, transform=ccrs.PlateCarree())
ax[0,1].coastlines()
fig.colorbar(cm01, ax=ax[0,1])

cm10 = ax[1,0].contourf(x, y, HC2_regress, 100, cmap=cmap1, transform=ccrs.PlateCarree())
ax[1,0].coastlines()
cm10.set_clim(-1*np.max(HC2_regress), np.max(HC2_regress))
fig.colorbar(cm10, label='$W/m^2$', ax=ax[1,0])

cm11 = ax[1,1].contourf(x, y, cov_L2, 10, cmap=cmap2, transform=ccrs.PlateCarree())
ax[1,1].coastlines()
fig.colorbar(cm11, ax=ax[1,1])

cm20 = ax[2,0].contourf(x, y, HC3_regress, 100, cmap=cmap1, transform=ccrs.PlateCarree())
ax[2,0].coastlines()
cm20.set_clim(-1*np.max(HC3_regress), np.max(HC3_regress))
fig.colorbar(cm20, label='$W/m^2$', ax=ax[2,0])

cm21 = ax[2,1].contourf(x, y, cov_L3, 10, cmap=cmap2, transform=ccrs.PlateCarree())
ax[2,1].coastlines()
fig.colorbar(cm21, ax=ax[2,1])

plt.savefig('/rds/general/user/jc9017/home/heatmapsNEW.png')

#%%
#Plot 7: Time Series per Basin

## Masking Heat Content Values ###
HC_ = np.where(np.isnan(HC), None, HC)
HC__ = ma.masked_where(HC_ == None, HC_)
## Indices ##
north_atl_X, north_atl_Y = [96, 185], [83, 152]
south_atl_X, south_atl_Y = [110, 199], [23,83]
north_pac_X, north_pac_Y = [0,99], [83,141]
south_pac_X, south_pac_Y = [0,99], [23,83]
ind_X, ind_Y = [199,299], [23,93]
antarc_X, antarc_Y = [0,360], [0,23]
arc_X, arc_Y = [0,360], [152,172]

## North Atlantic ##
HC_natl = []
for t in range(828):
    temp = []
    A=0
    for i in range(north_atl_X[0], north_atl_X[1]):
        for j in range(north_atl_Y[0], north_atl_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_natl.append(np.sum(temp)/A)
HC_natl = np.array(HC_natl) - np.mean(HC_natl)
natl_decomp = sm.tsa.seasonal_decompose(HC_natl, model='additive', freq = 12)
## South Atlantic ##
HC_satl = []
for t in range(828):
    temp = []
    A=0
    for i in range(south_atl_X[0], south_atl_X[1]):
        for j in range(south_atl_Y[0], south_atl_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_satl.append(np.sum(temp)/A)
HC_satl = np.array(HC_satl) - np.mean(HC_satl)
satl_decomp = sm.tsa.seasonal_decompose(HC_satl, model='additive', freq = 12)
## North Pacific ##
HC_npac = []
for t in range(828):
    temp = []
    A=0
    for i in range(north_pac_X[0], north_pac_X[1]):
        for j in range(north_pac_Y[0], north_pac_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_npac.append(np.sum(temp)/A)
HC_npac = np.array(HC_npac) - np.mean(HC_npac)
npac_decomp = sm.tsa.seasonal_decompose(HC_npac, model='additive', freq = 12)
## South Pacific ##
HC_spac = []
for t in range(828):
    temp = []
    A=0
    for i in range(south_pac_X[0], south_pac_X[1]):
        for j in range(south_pac_Y[0], south_pac_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_spac.append(np.sum(temp)/A)
HC_spac = np.array(HC_spac) - np.mean(HC_spac)
spac_decomp = sm.tsa.seasonal_decompose(HC_spac, model='additive', freq = 12)
## Indian Ocean ##
HC_ind = []
for t in range(828):
    temp = []
    A=0
    for i in range(ind_X[0], ind_X[1]):
        for j in range(ind_Y[0], ind_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_ind.append(np.sum(temp)/A)
HC_ind = np.array(HC_ind) - np.mean(HC_ind)
ind_decomp = sm.tsa.seasonal_decompose(HC_ind, model='additive', freq = 12)
## Arctic Ocean ##
HC_arc = []
for t in range(828):
    temp = []
    A=0
    for i in range(arc_X[0], arc_X[1]):
        for j in range(arc_Y[0], arc_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_arc.append(np.sum(temp)/A)
HC_arc = np.array(HC_arc) - np.mean(HC_arc)
arc_decomp = sm.tsa.seasonal_decompose(HC_arc, model='additive', freq = 12)
## Antarctic Ocean ##
HC_antarc = []
for t in range(828):
    temp = []
    A=0
    for i in range(antarc_X[0], antarc_X[1]):
        for j in range(antarc_Y[0], antarc_Y[1]):
            if not ma.is_masked(HC__[j,i,t]):
                j_ = y[j]*((np.pi)/180)
                A += ((np.pi/180)**2)*((6730e3)**2)*np.cos(j_)
                temp.append(HC__[j,i,t])
    HC_antarc.append(np.sum(temp)/A)
HC_antarc = np.array(HC_antarc) - np.mean(HC_antarc)
antarc_decomp = sm.tsa.seasonal_decompose(HC_antarc, model='additive', freq = 12)

#%%
plt.clf()
time = np.linspace(1,828,828)
fig, (ax1, ax2) = plt.subplots(2,1)
fig.set_size_inches(10,5)
ax1.plot(time_, natl_decomp.trend, label='North Atlantic')
ax2.plot(time_, satl_decomp.trend, label='South Atlantic')
ax1.plot(time_, npac_decomp.trend, label='North Pacific')
ax2.plot(time_, spac_decomp.trend, label='South Pacific')
ax2.plot(time_, ind_decomp.trend, label='Indian Ocean')
ax1.plot(time_, arc_decomp.trend, label='Arctic Ocean')
ax2.plot(time_, antarc_decomp.trend, label='Antarctic Ocean')
ax1.legend()
ax1.set_ylabel('Heat Content per unit Area ($J/{m^2}$)')
ax1.set_xlabel('Time in Years')
ax2.legend()
ax2.set_ylabel('Heat Content per unit Area ($J/{m^2}$)')
ax2.set_xlabel('Time in Years')
plt.savefig('/rds/general/user/jc9017/home/Basin_Efficiency_AnomalyNEW.png')

#%%

# PLot 8: Depth Penetration per Basin

def linear_fit(x,a,b):
    y = a*x + b
    return y
T_p0 = [300,250]
time = np.linspace(1,828,828)

natl_trends = []
for p in range(42):
    trends = []
    for i in range(north_atl_X[0], north_atl_X[1]):
        for j in range(north_atl_Y[0], north_atl_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    natl_trends.append(np.mean(trends))

satl_trends = []
for p in range(42):
    trends = []
    for i in range(south_atl_X[0], south_atl_X[1]):
        for j in range(south_atl_Y[0], south_atl_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    satl_trends.append(np.mean(trends))
    
npac_trends = []
for p in range(42):
    trends = []
    for i in range(north_pac_X[0], north_pac_X[1]):
        for j in range(north_pac_Y[0], north_pac_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    npac_trends.append(np.mean(trends))

spac_trends = []
for p in range(42):
    trends = []
    for i in range(south_pac_X[0], south_pac_X[1]):
        for j in range(south_pac_Y[0], south_pac_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    spac_trends.append(np.mean(trends))
    
ind_trends = []
for p in range(42):
    trends = []
    for i in range(ind_X[0], ind_X[1]):
        for j in range(ind_Y[0], ind_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    ind_trends.append(np.mean(trends))
    
arc_trends = []
for p in range(42):
    trends = []
    for i in range(arc_X[0], arc_X[1]):
        for j in range(arc_Y[0], arc_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    arc_trends.append(np.mean(trends))
    
antarc_trends = []
for p in range(42):
    trends = []
    for i in range(antarc_X[0], antarc_X[1]):
        for j in range(antarc_Y[0], antarc_Y[1]):
            if not ma.is_masked(thetaF[j,i,p,0]):
                T_param, T_cov = curve_fit(linear_fit, time, thetaF[j,i,p,:].astype(np.float64), p0 = T_p0)
                trends.append(T_param[0])
                
    antarc_trends.append(np.mean(trends))
            
#%%
for i in range(42):
    natl_trends[i] = 350*natl_trends[i] 
    satl_trends[i] = 350*satl_trends[i] 
    npac_trends[i] = 350*npac_trends[i] 
    spac_trends[i] = 350*spac_trends[i] 
    ind_trends[i] = 350*ind_trends[i] 
    arc_trends[i] = 350*arc_trends[i] 
    antarc_trends[i] = 350*antarc_trends[i] 
#%%
plt.clf()
fig, ax = plt.subplots(nrows=3, ncols=3, sharex=False, sharey=False)
fig.set_size_inches(14,10)
ax[0,0].plot(natl_trends, -z, '-', color='aquamarine')
ax[0,0].plot(natl_trends, -z, 'o', ms=2, color='cornflowerblue')
ax[0,0].text(0.6,0.1,'North Atlantic', transform=ax[0,0].transAxes)
ax[0,0].axvline(x=0, color='k', alpha=0.2)

ax[0,1].plot(npac_trends, -z, '-', color='aquamarine')
ax[0,1].plot(npac_trends, -z, 'o',ms=2, color='cornflowerblue')
ax[0,1].text(0.6,0.1,'North Pacific', transform=ax[0,1].transAxes)
ax[0,1].axvline(x=0, color='k', alpha=0.2)

ax[0,2].plot(arc_trends, -z, '-', color='aquamarine')
ax[0,2].plot(arc_trends, -z, 'o', ms=2, color='cornflowerblue')
ax[0,2].text(0.6,0.1,'Arctic Ocean', transform=ax[0,2].transAxes)
ax[0,2].axvline(x=0, color='k', alpha=0.2)

ax[1,1].plot(ind_trends, -z, '-', color='limegreen')
ax[1,1].plot(ind_trends, -z, 'o', ms=2, color='forestgreen')
ax[1,1].text(0.6,0.1,'Indian Ocean', transform=ax[1,1].transAxes)
ax[1,1].axvline(x=0, color='k', alpha=0.2)

ax[2,0].plot(satl_trends, -z, '-', color='salmon')
ax[2,0].plot(satl_trends, -z, 'o', ms=2, color='crimson')
ax[2,0].text(0.6,0.1,'South Atlantic', transform=ax[2,0].transAxes)
ax[2,0].axvline(x=0, color='k', alpha=0.2)

ax[2,1].plot(spac_trends, -z, '-', color='salmon')
ax[2,1].plot(spac_trends, -z, 'o', ms=2, color='crimson')
ax[2,1].text(0.6,0.1,'South Pacific', transform=ax[2,1].transAxes)
ax[2,1].axvline(x=0, color='k', alpha=0.2)

ax[2,2].plot(antarc_trends, -z, '-', color='salmon')
ax[2,2].plot(antarc_trends, -z, 'o', ms=2, color='crimson')
ax[2,2].text(0.1,0.1,'Antarctic Ocean', transform=ax[2,2].transAxes)
ax[2,2].axvline(x=0, color='k', alpha=0.2)

fig.text(0.5, 0.04, fontsize=18,s='Temperature Change per Annum (K/yr)', ha='center')
fig.text(0.08, 0.5, fontsize=18,s='Depth (m)', va='center', rotation='vertical')
fig.text(0.5, 0.9, fontsize=20, s='Plots of Depth Varying Temperature Change for Different Basins', ha='center')

plt.savefig('/rds/general/user/jc9017/home/Basin_Heat_PenetrationNEW.png')

#%%

            














