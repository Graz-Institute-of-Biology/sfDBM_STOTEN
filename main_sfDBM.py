# -*- coding: utf-8 -*-
"""
The simplieifed Field-scale desert biocrust model (sfDBM)

Reference:
Unravelling the main mechanism responsible for nocturnal CO2 uptake by dryland soils 
Kim et al. 
Science of the Total Environment
DOI:    

This Script simulates the soil-atmosphere CO2 exchange using envrionmental input (soil properties, temperatures, atmospheric CO2 etc.)

Additional CO2 reaction terms can be added with the unit of XXXX 

@author: Minsu Kim (minsu.kim@uni-graz.at) https://orcid.org/0000-0002-3942-3743

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from scipy.integrate import solve_ivp
import os
from functions import * # load functions used in this script

# %%
def main_tabernas_biocrusts(PATH, crust_name, realTinput, respif, fname_target, filesave=True, plotFigures=True,figuresave=False):
#
    df_soil_r, df_obs, averageT = load_field_data_soil_properties(PATH, crust_name, realTinput)
    
    porosity = df_soil_r.porosity.values[0]
    amountCations = 0.003+df_soil_r.Z_in.values[0] # Adjust Z values to match pH values
    initSulfate= 960 #960 mg/L ~ 10mM from MC crusts averaged values 
    initCalcium= 400 # 400mg/L ~ 10mM from MC crusts averaged values
    
    # % Set initial conditions of substrates, the steady state calculation based on the pH values at the saturated condition
    initC = steady_state(df_soil_r, averageT, porosity, amountCations,initSulfate,initCalcium)
    Zion = amountCations
    # mean soil weight per patch --- will be used to connect calcite and gypsum amount
    averageMsoil = 1
    # initial amount of solid calcite contents as a fraction
    initCalcitef = df_soil_r['CaCO3 (g/g)'].values[0] 
    # give input as gypsum contents
    initGypsumf = df_soil_r['gypsum (%)'].values[0]*0.01 
    # unit in g for solid phase on the surface
    initCalcite = averageMsoil*initCalcitef
    # unit in g for solid phase on the surface
    initGypsum = averageMsoil*initGypsumf
    # This value is the scaler for the CO2 dissolution affected by other substances (e.g. enzymes)
    # However, this work does not include this factor, thus set it to 1.
    fCA = 1

    dt = 60 #second %1 min.
    plottt = 10 # 10 minutes 
    timeLine = df_obs.index.to_list()
    examineH = timeLine[-1]
    dataPoints = len(timeLine)
    DeltaT = dt

    soilChemConditions = dict(
                              amountCations=amountCations,
                              Zion=Zion,
                              initCalcium=initCalcium,
                              initCalcite=initCalcite,
                              initGypsum=initGypsum,
                              initSulfate=initSulfate
                              )
    
    list_compounds = load_list_compounds_gypsum(soilChemConditions, averageT)
    
    MumaxTt,HenryConstList,pKList,Density_air = EnvironmentProfileNItrification(list_compounds, np.array(averageT),C_only=False)
    pKs = pKList.mean(axis=0).mean(axis=0)

    nx = 1
    ny = 1

    patchArea = df_soil_r['total_soil_SSA (m2/g)'].values[0] # in m^2/g -> The value indicates the ss area of the domain, thus [m2] 
    porosity = df_soil_r.porosity.values[0]
    waterFilm = df_obs['wft'].iloc[0]
    # Calculation of the unit volume of 1g of soil
    # Bulk density assumed to be 1.8 g/cm3 
    # total volume of the domain V = 1/1.8  = 0.55555 for a cube (1/1.8)**(1/3)=0.822cm length
    # From the porosity, calculate the pore volume and soil volume 
    poreVolume = porosity*(1/1.8)*10**-6 # in m^3 
    #Calculate the steady state for half saturation 
    waterVolume = waterFilm*patchArea # in m^3 
    gasVolume = poreVolume-waterVolume # in m^3 
    #Based on the total soil specific area, get the effective waterfilm thickness 
    voidTh = poreVolume/patchArea#Effective void thickness 
    
    #Reactive surface area of calicte and gypsum
    SA_calcite = initCalcite*df_soil_r['SSA (m2 g-1)'].values[0] #m2 
    SA_gypsum = patchArea*initGypsumf #m2
    
    ions = np.array(list_compounds.e_val.values)
    chemlist = list_compounds.index.values
    chemlist = np.delete(chemlist,[16])
    ions = np.delete(ions,[16])
    chemlist_woH = np.delete(chemlist,[15])

    # update initial conditions of concentration 
    for i, chem in enumerate(list_compounds.index):
        list_compounds.loc[chem,'initC'] = initC[chem].item()
        
    pH_ini = initC['pH'].item()

    # assign concentration and find mass equiiliria of gaseous substrates per patch according to Henry's law
    sitesCgas = dict()
    sitesNEquil = dict()
    sitesNGEquil = dict()
    HenryDomain = dict()
    effluxList = dict()
    for c_name in list(HenryConstList):
        sitesCgas.update({c_name:list_compounds.loc[c_name].initCg*np.ones((nx,ny))})
        sitesNGEquil.update({c_name:list_compounds.loc[c_name].partial_pressure*Density_air})
        HenryDomain.update({c_name:HenryConstList[c_name]})
        sitesNEquil.update({c_name:np.multiply(HenryConstList[c_name],sitesCgas[c_name])})
        effluxList.update({c_name:np.zeros((int(dataPoints),1))})
        
    sitesC = dict()
    sitesN = dict()
    nutShare = dict()
    totConsumed = dict()
    changeN = dict()
    sitesCtemp = dict()  # temporary arrays for pH calculations
    timeConcDist = dict()  # temporary arrays for pH calculations
    for chem in list_compounds.index:
        sitesC.update({chem: list_compounds.loc[chem].initC*np.ones((nx, ny))})
        sitesN.update(
            {chem: list_compounds.loc[chem].initC*waterVolume*np.ones((nx, ny))})
        nutShare.update({chem: np.ones((nx, ny))})
        totConsumed.update({chem: np.zeros((nx, ny))})
        changeN.update({chem: np.zeros((nx, ny))})
        sitesCtemp.update({chem: np.zeros((nx, ny))})
        timeConcDist.update({chem: np.zeros((nx, ny, dataPoints))})
    sitesC.update({'pH': pH_ini*np.ones((nx, ny))})
    timeConcDist.update({'pH': pH_ini*np.ones((nx, ny, dataPoints))})
    temperaturelist = np.zeros(dataPoints)

    # % Starting dynamics
    print('Day 0 starts')
    for examineT in range(dataPoints): #Time for saving the values (interval of the data)
        for tt in range(plottt): #dynamics within the time unit

            t = examineT*plottt + tt
            exposedHours = t * dt /3600
            if exposedHours%6 == 0:
                print(exposedHours)   
            
            waterFilm = df_obs['wft'].iloc[examineT]
            waterVolume = patchArea*waterFilm
            gasVolume = (voidTh-waterFilm)*patchArea
        
            totalweight = waterVolume*998000*np.ones((nx, ny))  # weight of water in g
            for chem in list_compounds.index:
                sitesC[chem] = sitesN[chem]/waterVolume
                if not chem in ['CaCO3', 'CaSO42H2O']:
                    # weight of soil solution
                    totalweight = totalweight + sitesN[chem].item()
        

            thetaTC = df_obs['ts'].iloc[examineT]
                
            temperaturelist[examineT] = thetaTC
            # needs when temperature changes
            MumaxTt,HenryConstList,pKList,soil_air_density = EnvironmentProfileNItrification(list_compounds, np.array(thetaTC),C_only=False)
            Density_air = air_density(df_obs['ta'].iloc[examineT],0,101325.0) * 1000
            xc_atm = df_obs['xc_atm_pd'].iloc[examineT]/1000
            #update equilibrium conditions based on changes in temperature 
            for chem in list(HenryConstList):
                if chem == 'CO2':
                    sitesNGEquil.update({chem:xc_atm})                
                else:
                    sitesNGEquil.update({chem:list_compounds.loc[chem].partial_pressure*Density_air})
                HenryDomain.update({chem:HenryConstList[chem]})
                sitesNEquil.update({chem:np.multiply(HenryConstList[chem],sitesCgas[chem])})

            hlist = []
            gclist = []
            for chem in list(HenryConstList):
                hlist.append(HenryDomain[chem])
                gclist.append(sitesNGEquil[chem]*0.001/list_compounds.loc[chem].molwt)

            for xx in range(nx):
                for yy in range(ny):
                    
                    InvasedIsland =1
                    kla = 1.e-9*patchArea/waterVolume
                    gasflow = InvasedIsland*0.01/60
                    pKs = pKList[xx, yy, :]
                    theta = thetaTC
                    Z = sitesN['Zion'][xx, yy]/totalweight[xx, yy]
                    totalw = totalweight[xx, yy]
                    rlist = []
                    for chem in list(chemlist):

                        if chem in ['CO2']:
                            rlist.append(respif*0.1*5.e-6)
                        else:
                            rlist.append(0)
                    
                    # find equilibrium of the chemical conditions and update values.
                    y0 = []
                    for chem in list(chemlist):  # compounds not related to gypsum
                        if not chem in ['CaCO3', 'CaSO42H2O']:
                            y0.append(sitesN[chem][xx, yy]*1000/list_compounds.loc[chem].molwt/totalw)
                        else:
                            y0.append(sitesC[chem][xx, yy]/list_compounds.loc[chem].molwt/totalw) # in mol/g
                            
                    # y0[:len(ions)] : aqueuous dissolved substances
                    # add gas concentrations
                    for chem in list(HenryConstList):
                        y0.append(sitesCgas[chem][xx, yy]*0.001/list_compounds.loc[chem].molwt)
                     # y0[len(ions):] : gaseous substances
                    for chem in list(HenryConstList):  # efflux tracking
                        y0.append(0)

                    y0 = np.array(y0)
                    y0[y0 < 0] = 0
                    soln = solve_ivp(chem_gypsum_calcite_aw_Z_reaction, (0.0, DeltaT), y0,
                                     args=(Z, ions, pKs, hlist, gclist, rlist, kla, gasflow /
                                           gasVolume, theta, totalw, SA_gypsum, SA_calcite, fCA),
                                     method='BDF', rtol=1E-8, atol=1E-8)

                    solY = soln.y[:25, -1]
                    solY[solY < 0] = 0

                    for i, chem in enumerate(chemlist):
                        if not chem in ['CaCO3', 'CaSO42H2O']: #Precipitated calcite and gypsum
                            sitesCtemp[chem][xx, yy] = soln.y[i, -1]*0.001*list_compounds.loc[chem].molwt*totalw/waterVolume
                        else:
                            sitesCtemp[chem][xx, yy] = soln.y[i, -1]*list_compounds.loc[chem].molwt*totalw #back to g

                    for i, chem in enumerate(list(HenryConstList), start=len(ions)):
                        sitesCgas[chem][xx, yy] = solY[i]*1000*list_compounds.loc[chem].molwt

                    for i, chem in enumerate(list(HenryConstList), start=len(ions)+len(HenryConstList)):
                        effluxList[chem][examineT, 0] = effluxList[chem][examineT,0] + soln.y[i, -1]*1000*list_compounds.loc[chem].molwt*waterVolume # unit in g
    
            for chem in chemlist_woH:
                sitesC[chem] = sitesCtemp[chem]
            
            sitesC['H+'] = sitesCtemp['H+']
            sitesC['pH'] = -np.log10(sitesCtemp['H+']*0.001)
              
                
            for chem in list(sitesC)[:-1]:
                sitesN[chem] = sitesC[chem]*waterVolume

        
        for chem in list(sitesC)[:-1]:
            timeConcDist[chem][:,:,examineT] = sitesC[chem]
            sitesN[chem] = sitesC[chem]*waterVolume
            
        timeConcDist['pH'][:,:,examineT] = sitesC['pH']


    if plotFigures:    
        plot_results(fname_target,df_obs,effluxList,timeConcDist,timeLine,examineH,18,7)
          
    if filesave:
        dbm_dynamics = os.path.join(PATH, "model_output", fname_target+".npz")
        np.savez(dbm_dynamics,
                  examineH=examineH,
                  timeLine=timeLine,
                  chemlist=chemlist,
                  timeConcDist=timeConcDist,
                  effluxList=effluxList,
                  soilChemConditions=soilChemConditions,
                  temperaturelist = temperaturelist,
                  )        

    return timeLine,effluxList, timeConcDist

# %%

PATH = "C:/Users/kimmi/Documents/GitHub/dfDBM_STOTEN"

filesave = True
figuresave = False
plotFigures = True

crust_name = 'MC' # ['PD', 'MC', 'SD']
realTinput = 'winter' #['winter', 'summer']
respif = 0.0

fname_target = 'sfDBM_'+crust_name+'_'+realTinput+'_with_reac_respif_%.4f'%respif 
print(fname_target)

try:
    timeLine,effluxList, timeConcDist = main_tabernas_biocrusts(PATH, crust_name, realTinput,respif, fname_target, filesave=filesave, plotFigures=plotFigures,figuresave=figuresave)
except:
    print('error occured')

#%%
#if __name__ == '__main__':
#    main()