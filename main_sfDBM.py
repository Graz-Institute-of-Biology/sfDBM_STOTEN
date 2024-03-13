# -*- coding: utf-8 -*-
"""
The simplified Field-scale desert biocrust model (sfDBM)

Reference:
Unravelling the main mechanism responsible for nocturnal CO2 uptake by dryland soils 
Kim et al.(2024) 
Science of the Total Environment
DOI:    

This Script simulates the CO2 efflux from biocrust soils using previously observed envrionmental input (soil properties, temperatures, atmospheric CO2 etc.)
Additional CO2 reaction terms can be added with the unit of XXXX 

@author: Minsu Kim (minsu.kim@uni-graz.at) https://orcid.org/0000-0002-3942-3743

"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob
from scipy.integrate import solve_ivp
import os
import argparse
import textwrap

# %%

def crust_to_repi(crust_name):
    """
    The function to call the example input data (PD, MC, and SD crusts)
    """   
    switcher = {
        'PD': 1,
        'MC': 7,
        'SD': 12,
    }
    return switcher.get(crust_name, "no data provided ")

def load_field_data_soil_properties(PATH, crust_name, realTinput):
    """
    The function to load soil properties and field observation data
    """    
    
    # load pre-calculated water retention and water film thickness (Kim and Or 2016)
    # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0147394
    preRetCV = os.path.join(PATH, "input_data", "Retention_hydration_interp_pot_D_2.65.csv")
    df_ret = pd.read_csv(preRetCV)
    df_ret['local_satu'] = df_ret.FilmThickTemp/df_ret.FilmThickTemp.iloc[0]

    if realTinput == 'summer':
        t1 = '2019-07-01'
        t2 = '2019-08-01'
        avgT = 35
        
    elif realTinput == 'winter':
        t1 = '2019-01-01'
        t2 = '2019-02-01'
        avgT = 10
        
    else:
        t1 = ''
        t2 = ''
        avgT = 25
        
    rep_i = crust_to_repi(crust_name)
    
    #Load soil properties and pre-calulated input of Zin 
    soilprop = os.path.join(PATH, "input_data", "data_replicates"+"_avgT_%d"%avgT+"_zin.csv")
    df_soil = pd.read_csv(soilprop)
    df_soil_r = df_soil[df_soil.replicate == rep_i]

    # Load the example time sereies data observed at the field (winter and summer time data are provided)
    # Details fo the field data can be found at (Lopez-Canfin et al 2022)
    # https://www.sciencedirect.com/science/article/pii/S0016706122003743
    fname = 'time_series_*_'+t1+'_'+t2+'_rep_%d'%int(rep_i)
    fnames = glob.glob(os.path.join(PATH, "input_data", fname+".csv"))
    tempList = pd.read_csv(fnames[0])
            
    tempList = tempList.set_index(pd.to_datetime(tempList.time))        
    tempList['intp_time'] = (tempList.index - tempList.index[0]).total_seconds()/3600
    tempList = tempList.reset_index(drop=True).set_index('intp_time')
    timeL = np.arange(0,tempList.index.max()+1/6,1/6) 
    tempList = tempList.reindex(timeL).interpolate()# Change hourly data to 10mins data
    tempList = tempList.dropna()        
    porosity = df_soil_r['porosity'].values[0]
    
    # Based on the water content at the field condition, back calculate water film thickness 
    # water content and porosity -> saturation -> film thickness
    wfseries = []
    for swc in tempList.swc_md:
        ids = len(df_ret[df_ret.local_satu*porosity>swc])
        wfseries.append(np.interp(swc, [df_ret.iloc[ids].local_satu*porosity,df_ret.iloc[ids-1].local_satu*porosity],[df_ret.iloc[ids].FilmThickTemp,df_ret.iloc[ids-1].FilmThickTemp]))
    
    tempList['wft'] = wfseries

    return df_soil_r, tempList, avgT #Return soil property, field observation, average tempeature assigned for approximated chemical conditions


def load_list_compounds_gypsum(soilChemConditions, averageT):
    """
    The function to load parameters of chemical conditions (substrate concentrations, gas solubility, etc.)
    We note here that, this function is a sub-module of a general soil model for C and N cycle.
    Therefore, this function includes inorganic C and N, also the carbonates and gypsum, gaseous compounds relavant for the N cycle
    In the paper published at STOTEN (), we have focused on the DIC fractionation with an assumption that other compounds given here are assumed to be absent.
    This function can be modified and used in the other contexts
    """  
    #list of substrates 
    list_substrates = ['O2', # Oxygen 
                       'CO2',# Carbon dioxide
                       'HCO3-',# Bicarbonate
                       'CH2O', # Formaldehyde (assumed to be the simplist organic C compound)
                       'NO3-', # Nitrate
                       'NH3', # Ammonia 
                       'NO2-', # Nitrite
                       'HONO', # Nitrous acid
                       'NO', # Nitric oxide
                       'N2O', # Nitrous oxide
                       'NH4+', # Ammonium
                       'CO3 2-', #Carbonate
                       'Ca2+', # Calcium
                       'CaCO3*',#Complexed Ca2+-CO32-
                       'CaCO3', #Calcium carbonate
                       'H+', # Hydrogen (used for pH)
                       'Zion', # Assumed to be auxilary ions for charge balance 
                       'CaSO42H2O', #Gypsum
                       'SO42-', #Sulphate 
                       'CaSO4*' #Complexed Ca2+-SO42-
                       ] 
    #molecular weights of substrates (g/mol)
    list_Molwt = [32.00,44.01,61.02,30.03,62.00,17.03,46.01,47.013, 30.01,44.01,18.04,60.01,40.08,100.09,100.09,1.00, 1, 172.17, 96.07, 136.14] 
    
    # e valancy
    list_evalancy = [0,0,-1,0,-1,0,-1,0,0,0,1,-2,2,0,0,1,1,0,-2,0]

    # We note here that Zion indicates the amount of cations, 
    # which is similar to Ca2+ but does not react with inorganic carbon. 
    # The moelcular weight is assumed to be the same of Ca and valancy 1

    # Diffusion coefficients
    list_diffusion = np.array([1.73, 1.65, 1.02, 0.5184, 1.47, 1.41, 1.64, 1.64, 1.91, 1.58, 1.41, 0.79, 0.8640, 0, 0, 8.04, 0, 0, 0.92448, 0])*0.0001/(24*3600) #m^2/s
    # NO diffusivity from Zacharia Deen (2004) https://link.springer.com/content/pdf/10.1007/s10439-005-8980-9.pdf    
    # SO4 2- diffusivity from https://www.aqion.de/site/diffusion-coefficients
    # CaSO4 diffusivity assumed zero
    
    d = {'molwt':list_Molwt, 'e_val':list_evalancy, 'diffusivity':list_diffusion}
    list_compounds = pd.DataFrame(index=list_substrates, data=d)
    
    # List of gaseous compounds 
    list_gas_substances = ['O2', # Oxygen 
                           'CO2',# Carbon dioxide
                           'NH3', # Ammonia 
                           'HONO', # Nitrous acid
                           'NO', # Nitric oxide
                           'N2O' # Nitrous oxide
                           ]
    

    NH3ppb = 1
    HONOppb = 0
    NOppb = 0
    N2Oppb = 333
    CO2ppm = 400
    
    #partial pressure of atmospheric substances
    pN2 = 0.7809 #when N2 is considered for finxation
    pO2 = 0.2095
    pCO2 = CO2ppm * 10 ** (- 6) #assuemd to be 400 ppm
    pNH3 = NH3ppb * 10 ** (- 9)
    pHONO = HONOppb * 10 ** (- 9)
    pN2O = N2Oppb * 10 ** (- 9)
    pNO = NOppb * 10 ** (- 9)
    PartialPressure = [pO2,pCO2,pNH3,pHONO,pNO,pN2O]
    
    # parameters needed for calculating temperature-dependent Henry's solubility
    # All information from https://www.henrys-law.org/henry/ 
    # conversion to M/atm (1 mol/m3Pa = 101.325M/atm) 
    convH = 101.325
    #Oxygen KH Warneck and Williams (2012)
    KHO2 = convH*1.2 * 10 ** (- 5)
    enthalpyO2 = 1700
    #Carbondioxide KH  	Sander et al. (2011)
    KHCO2 = convH*3.3 * 10 ** (- 4)
    enthalpyCO2 = 2400
    #Ammonia KH  	Sander et al. (2011)
    KHNH3 = convH*5.9 * 10 ** (- 1)
    enthalpyNH3 = 4200
    #HONO KH Park and Lee (1988)
    KHHONO = convH*4.8 * 10 ** (- 1)
    enthalpyHONO = 4900
    #NO KH Warneck and Williams (2012)
    KHNO = convH*1.9 * 10 ** (- 5)
    enthalpyNO = 1600
    #N2O KH Warneck and Williams (2012)
    KHN2O = convH*2.4 * 10 ** (- 4)
    enthalpyN2O = 2700
    
    #Arrhenius based HONO dissociation
    KaHONO = 5.6 * 10 ** (-4);
    enthalpyHONOKa = 1423.5;
    
    list_gas_KH = [KHO2, KHCO2,KHNH3, KHHONO,KHNO, KHN2O]
    list_gas_enthalpy = [enthalpyO2, enthalpyCO2,enthalpyNH3, enthalpyHONO,enthalpyNO, enthalpyN2O]
    
    d = {'partial_pressure':PartialPressure, 'KH':list_gas_KH, 'enthalpy':list_gas_enthalpy}
    list_gas_compounds =pd.DataFrame(index=list_gas_substances, data=d)
    
    list_compounds = pd.concat([list_compounds, list_gas_compounds], axis=1)
    
    # Initial conditions given based on the soil chemical properies and observation
    amountCations = soilChemConditions['amountCations']
    Zion = soilChemConditions['Zion']
    initCalcium = soilChemConditions['initCalcium']
    initCalcite = soilChemConditions['initCalcite']
    initGypsum = soilChemConditions['initGypsum']
    initSulfate = soilChemConditions['initSulfate']

    #Dissociation constants of DIC at the stardard condition
    KA = 10 ** (- 9.4003)
    K1C = 10 ** (- 6.3819)
    K2C = 10 ** (- 10.3767)
    
    #Atmospheric pressure
    pa = 101325.0 # [Pa]
    rho = air_density(averageT,0,pa) * 1000 #air density in the unit of [kg/m^3] at the given T and Pa

    list_compounds['initC'] = 0
    list_compounds['initCg'] = 0

    list_compounds.loc['O2','initCg'] = rho*list_compounds.loc['O2'].partial_pressure #mg/L
    list_compounds.loc['CO2','initCg'] = rho*list_compounds.loc['CO2'].partial_pressure #mg/L 
    list_compounds.loc['NH3','initCg'] = rho*list_compounds.loc['NH3'].partial_pressure #mg/L  a few ppb: assumed 5 %Gong et al 2011: atmospheric level of ammonia
    list_compounds.loc['HONO','initCg'] = rho*list_compounds.loc['HONO'].partial_pressure #mg/L %1ppb for HONO
    list_compounds.loc['NO','initCg'] = rho*list_compounds.loc['NO'].partial_pressure #mg/L  a few ppb
    list_compounds.loc['N2O','initCg']= rho*5*list_compounds.loc['N2O'].partial_pressure #mg/L %500ppb (o.5ppm) for N2O
    list_compounds.loc['O2','initC'] = list_compounds.loc['O2','initCg']*solubility(list_compounds.loc['O2'].KH, list_compounds.loc['O2'].enthalpy, averageT) # N4 : O2 mg/L : Henry's law at 1 atm
    list_compounds.loc['CO2','initC'] = list_compounds.loc['CO2','initCg']*solubility(list_compounds.loc['CO2'].KH, list_compounds.loc['CO2'].enthalpy, averageT) # N1 : CO2 mg/L : Henry's law at 1 atm
    list_compounds.loc['HCO3-','initC'] = list_compounds.loc['CO2','initC']*61.0171/44.0098*K1C/10**(-7.5) # at equilbrium with atmospheric level
    list_compounds.loc['CH2O','initC'] = 26.9958*0.1 # N5 : SS 0.1M (readily degradable organic substrate) : CH1.5O0.5N0.1 %Initmole only controls the sugar concentraion. (Catillo-Monroy 2010): 34.61mg/kg soil: with porosity 0.3: 34.61*2.6*0.3~
    list_compounds.loc['NO3-','initC'] = 0.6 #assumed value, does not play a role for abiotic partinioning of DIC
    list_compounds.loc['NO2-','initC'] = 0 #assumed to be zero
    list_compounds.loc['HONO','initC'] = list_compounds.loc['HONO','initCg']*solubility(list_compounds.loc['HONO'].KH, list_compounds.loc['HONO'].enthalpy, averageT) #HONO in solution
    list_compounds.loc['NO','initC'] = 0  #NO in solution assumed to be zero
    list_compounds.loc['N2O','initC'] = list_compounds.loc['N2O','initCg']*solubility(list_compounds.loc['N2O'].KH, list_compounds.loc['N2O'].enthalpy, averageT) #N2O
    list_compounds.loc['NH3','initC'] = list_compounds.loc['NH3','initCg']*solubility(list_compounds.loc['NH3'].KH, list_compounds.loc['NH3'].enthalpy, averageT) #NH3 (almost zero)
    list_compounds.loc['Ca2+','initC'] = initCalcium #Ca2+ (calcium ion) in M (Johnson 2005 relevance)
    list_compounds.loc['H+','initC'] = 1000*10**(-7.5) #[H+] assumed to bepH 7.5 in mg/L
    list_compounds.loc['CO3 2-','initC'] = list_compounds.loc['HCO3-','initC']*60.010/61.0171*K2C/10**(-7.5) #CO3 2- based on the pH values 
    # calcite complex / precipitation
    list_compounds.loc['CaCO3*','initC'] = list_compounds.loc['Ca2+','initC']*list_compounds.loc['CO3 2-','initC']*10 ** (3.22)*0.001*100.090/(40.080*60.010) #Ca2+ (calcium ion) in M (Johnson 2005 relevance)
    list_compounds.loc['CaCO3','initC'] = initCalcite #amount of calciate in the unit of g ---> this will be transformed to mol in chemical module
    list_compounds.loc['NH4+','initC'] = 0
    list_compounds.loc['Zion','initC'] = Zion #[H+] in M
    list_compounds.loc['CaSO42H2O','initC'] = initGypsum #amount of gypsum on surface in the unit of g 
    list_compounds.loc['SO42-','initC'] = initSulfate #Sulfate ion in mg/L
    list_compounds.loc['CaSO4*','initC'] = initCalcium*initSulfate*4.9*10**(-6)*0.001*136.140/(40.080*96.070) #at equilibrium

    return list_compounds # Return list of chemical parameters and initial conditions for each compounds


def air_density(t = None,hr = None,p = None): 
    """
    Air density as a function of temperature and pressure 
    """  
# 1)'Equation for the Determination of the Density of Moist Air' P. Giacomo  Metrologia 18, 33-40 (1982)
# 2)'Equation for the Determination of the Density of Moist Air' R. S. Davis Metrologia 29, 67-70 (1992)
    T0 = 273.15
    T = T0 + t
    # 1) Coefficients values
    R = 8.31451   
    Mv = 18.015 * 10 ** - 3
    Ma = 28.9635 * 10 ** - 3
    A = 1.2378847 * 10 ** - 5
    B = - 1.9121316 * 10 ** - 2
    C = 33.93711047
    D = - 6.3431645 * 10 ** 3
    a0 = 1.58123 * 10 ** - 6
    a1 = - 2.9331 * 10 ** - 8
    a2 = 1.1043 * 10 ** - 10
    b0 = 5.707 * 10 ** - 6
    b1 = - 2.051 * 10 ** - 8
    c0 = 1.9898 * 10 ** - 4
    c1 = - 2.376 * 10 ** - 6
    d = 1.83 * 10 ** - 11
    e = - 0.765 * 10 ** - 8
    #-------------------------------------------------------------------------
# 2) Calculation of the saturation vapour pressure at ambient temperature, in [Pa]
    psv = np.exp(np.multiply(A,(T ** 2)) + np.multiply(B,T) + C + D / T)
    #-------------------------------------------------------------------------
# 3) Calculation of the enhancement factor at ambient temperature and pressure
    fpt = 1.00062 + (3.14 * 10 ** - 8) * p + (5.6 * 10 ** - 7) * (t ** 2)
    #-------------------------------------------------------------------------
# 4) Calculation of the mole fraction of water vapour
    xv = np.multiply(np.multiply(np.multiply(hr,fpt),psv),(1.0 / p)) * (10 ** - 2)
    #-------------------------------------------------------------------------
# 5) Calculation of the compressibility factor of air
    Z = 1 - (np.multiply((p / T),(a0 + a1 * t + a2 * (t ** 2) + np.multiply((b0 + b1 * t),xv) + np.multiply((c0 + c1 * t),(xv ** 2))))) + (np.multiply((p ** 2 / T ** 2),(d + np.multiply(e,(xv ** 2)))))
    #-------------------------------------------------------------------------
# 6) Final calculation of the air density in [kg/m^3]
    ro = np.multiply((np.multiply(p,Ma) / (np.multiply(np.multiply(Z,R),T))),(1 - np.multiply(xv,(1 - Mv / Ma))))
    return ro # [kg/m^3]


def solubility(KH, enthalpy, temperatureTemp):
    """
    Temperature dependent solubility
    """  
    return (temperatureTemp + 273.15) * KH * np.exp(enthalpy * ((1 / (temperatureTemp + 273.15)) - (1 /  298.15))) / 12.2
    


def steady_state(df_soil_r, averageT,porosity=0.5, amountCations=0.003,initSulfate=960,initCalcium= 400):
    """
    This function finds an initial condition for the field application. 
    Assuming the given temperature (no dynamics) and provied chemical properties, it calculates the pationing of substrates at the saturated condtions
    Here, the saturation is assuemd to be 0.99 for a numerical reason
    The function returns the initial concentration of the chemical substrates, that are in equlibilrium with gas-liquid partionioning at the given pH
    """      
    #This fCA factor is to accelerate the CO2 dissolution mimicking the role of carbonic anhydrase
    #Setting this value 1, indicates no effects in CO2 dissolution 
    fCA=1 

    # initial condition for chemical soil properites    
    initCalcitef = df_soil_r['CaCO3 (g/g)'].values[0]
    initGypsumf = df_soil_r['gypsum (%)'].values[0]*0.01
    target_pH = df_soil_r['pH'].values[0]
    Zion = amountCations
    # mean soil weight per patch --- will be used to connect calcite and gypsum amount
    # calculation for 1 g of soil
    averageMsoil = 1
    # unit in g for solid phase of the domain
    initCalcite = averageMsoil*initCalcitef
    # unit in g for solid phase of the domain
    initGypsum = averageMsoil*initGypsumf

    soilChemConditions = dict(
                              amountCations=amountCations,
                              Zion=Zion,
                              initCalcium=initCalcium,
                              initCalcite=initCalcite,
                              initGypsum=initGypsum,
                              initSulfate=initSulfate
                              )
    # Load all the compounds
    list_compounds = load_list_compounds_gypsum(soilChemConditions, averageT)
    
    MumaxTt, HenryConstList, pKList, Density_air = EnvironmentProfileNItrification(list_compounds, np.array(averageT), C_only=False)
    pKs = pKList.mean(axis=0).mean(axis=0)
    
    nx = 1
    ny = 1
    
    saturation = 0.99
    # porosity = df_soil_r.porosity.values[0]
    # Calculation of the unit volume of 1g of soil
    # Bulk density assumed to be 1.8 g/cm3 
    # total volume of the domain V = 1/1.8  = 0.55555 for a cube (1/1.8)**(1/3)=0.822cm length
    # From the porosity, calculate the pore volume and soil volume 
    # porevolume = \phi V = df_soil_r.porosity.values[0]*(1/1.8)
    poreVolume = porosity*(1/1.8)*10**-6 # in m^3 
    #Calculate the steady state for half saturation 
    waterVolume = saturation*poreVolume # in m^3 
    gasVolume = poreVolume-waterVolume # in m^3 
    #Based on the total soil specific area, get the effective waterfilm thickness 
    ssArea = df_soil_r['total_soil_SSA (m2/g)'].values[0] # in m^2/g -> The value indicates the ss area of the domain, thus [m2] 
    patchArea = ssArea
    waterFilm = waterVolume/ssArea # unit in [m]
    voidTh = poreVolume/ssArea#Effective void thickness 
    
    #Reactive surface for calcite and gypsum
    SA_calcite = initCalcite*df_soil_r['SSA (m2 g-1)'].values[0] #m2 
    SA_gypsum = ssArea*initGypsumf #m2
    
    # total weight of solution in g base calculation for an inital condition
    totalweight = waterVolume*998000
    # find equilibrium of the chemical conditions and update values.
    ions = np.array(list_compounds.e_val.values)
    chemlist = list_compounds.index.values
    chemlist = np.delete(chemlist,[16])
    ions = np.delete(ions,[16])
    chemlist_woH = np.delete(chemlist,[15])
    
    # # add gas concentrations
    hlist = []
    gclist = []
    for chem in list(HenryConstList):
        hlist.append(HenryConstList[chem])
        gclist.append(
            list_compounds.loc[chem].initCg*10**(-3)/list_compounds.loc[chem].molwt)
    
    pH_ini = target_pH
    nx = 1
    ny = 1
    
    # assign concentration and find mass equiiliria of gaseous substrates per patch according to Henry's law
    sitesCgas = dict()
    sitesNEquil = dict()
    sitesNGEquil = dict()
    HenryDomain = dict()
    AlphaM = dict()
    
    for c_name in list(HenryConstList):
        sitesCgas.update({c_name:list_compounds.loc[c_name].initCg*np.ones((nx,ny))})
        sitesNGEquil.update({c_name:list_compounds.loc[c_name].partial_pressure*Density_air})
        HenryDomain.update({c_name:HenryConstList[c_name]})
        sitesNEquil.update({c_name:np.multiply(HenryConstList[c_name],sitesCgas[c_name])})
        AlphaM.update({c_name:np.multiply(HenryConstList[c_name],waterVolume) + gasVolume})
    
    
    sitesC = dict()
    sitesN = dict()
    sitesCtemp = dict()  # temporary arrays for pH calculations
    
    for chem in list_compounds.index:
        sitesC.update({chem: list_compounds.loc[chem].initC*np.ones((nx, ny))})
        sitesN.update(
            {chem: list_compounds.loc[chem].initC*waterVolume*np.ones((nx, ny))})
        sitesCtemp.update({chem: np.zeros((nx, ny))})
    
    sitesC.update({'pH': pH_ini*np.ones((nx, ny))})
        
    #Here, the total weight indicates soil water weight not the soil weight
    totalweight = waterVolume*998000*np.ones((nx, ny))  # weight of water in g
    
    for chem in list_compounds.index:
        if not chem in ['CaCO3', 'CaSO42H2O']:
            totalweight = totalweight + sitesN[chem]
    
    # needs when temperature changes
    MumaxTt, HenryConstList, pKList, Density_air = EnvironmentProfileNItrification(
        list_compounds, np.array(averageT), C_only=False)
    
    # update equilibrium conditions based on changes in temperature
    for chem in list(HenryConstList):
        sitesNGEquil.update(
            {chem: list_compounds.loc[chem].partial_pressure*Density_air})
        HenryDomain.update({chem: HenryConstList[chem]})
    
    
    hlist = []
    gclist = []
    for chem in list(HenryConstList):
        hlist.append(HenryDomain[chem])
        gclist.append(sitesNGEquil[chem]*0.001/list_compounds.loc[chem].molwt)
    
    xx = 0
    yy = 0 
    
    kla = 1.e-9*patchArea/waterVolume
    gasflow = 100000
    
    pKs = pKList[xx, yy, :]
    theta = averageT
    totalw = totalweight[xx, yy]
    Z = Zion

    rlist = []
    for chem in list(chemlist):
        rlist.append(0)
    # find equilibrium of the chemical conditions and update values.
    y0 = []
    for chem in list(chemlist):  # compounds not related to gypsum
        #y0.append(sitesC[chem][xx,yy]*10**(-3)/list_compounds.loc[chem].molwt) # unit in M (mol/L)
        # unit in mol/kg
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

    soln = solve_ivp(chem_gypsum_calcite_aw_Z_reaction, (0.0, 1000), y0,
                     args=(Z, ions, pKs, hlist, gclist, rlist, kla, gasflow /
                           gasVolume, theta, totalw, SA_gypsum, SA_calcite, fCA),
                     method='BDF', rtol=1E-8, atol=1E-8)

    solY = soln.y[:25, -1]
    solY[solY < 0] = 0
    
    for i, chem in enumerate(chemlist):
        if not chem in ['CaCO3', 'CaSO42H2O']:
            sitesCtemp[chem][xx, yy] = soln.y[i, -1] *0.001*list_compounds.loc[chem].molwt*totalw/waterVolume
        else:
            sitesCtemp[chem][xx, yy] = soln.y[i, -1]*list_compounds.loc[chem].molwt*totalw #back to g
               
    for chem in chemlist_woH:
        sitesC[chem] = sitesCtemp[chem]
    
    sitesC['H+'] = sitesCtemp['H+']
    sitesC['pH'] = -np.log10(sitesCtemp['H+'])

    return sitesC

def EnvironmentProfileNItrification(list_compounds, Tdist,C_only): 
    
    stdTemp = 298.15
    absoluteT = 273.15
    conversionC = 12.2
    
    Tl = 297.7
    Th = 314.7
    Hh = 687900
    Hl = - 141100
    Ha = - 5430
    pa = 101325.0
    
    #HONO dissociation
    KaHONO = 5.6 * 10 ** (- 4)
    enthalpyHONOKa = 1423.5
    enthalpyHONOoverR = enthalpyHONOKa * 1000 / 8.314

    T = Tdist + absoluteT
    MumaxT = (1 / stdTemp) * T * np.exp((Ha / 8.31) * (1 / stdTemp - 1 / T)) / (1 + np.exp((Hl / 8.31) * (1 / Tl - (1 / T))) + np.exp((Hh / 8.31) * (1 / Th - (1 / T))))
    
    #Oxygen KH
    KHO2 = list_compounds.loc['O2'].KH
    enthalpyO2 = list_compounds.loc['O2'].enthalpy
    #Carbondioxide KH
    KHCO2 = list_compounds.loc['CO2'].KH
    enthalpyCO2 = list_compounds.loc['CO2'].enthalpy
    HenryConstListO2 = T * KHO2 * np.exp(enthalpyO2 * ((1 / T) - (1 / stdTemp))) / conversionC
    HenryConstListCO2 = T * KHCO2 * np.exp(enthalpyCO2 * ((1 / T) - (1 / stdTemp))) / conversionC

    if not C_only:
        #Ammonia KH
        KHNH3 = list_compounds.loc['NH3'].KH
        enthalpyNH3 = list_compounds.loc['NH3'].enthalpy
        #HONO KH
        KHHONO = list_compounds.loc['HONO'].KH
        enthalpyHONO = list_compounds.loc['HONO'].enthalpy
        #NO KH
        KHNO = list_compounds.loc['NO'].KH
        enthalpyNO = list_compounds.loc['NO'].enthalpy
        #N2O KH
        KHN2O = list_compounds.loc['N2O'].KH
        enthalpyN2O = list_compounds.loc['N2O'].enthalpy
        HenryConstListNH3 = T * KHNH3 * np.exp(enthalpyNH3 * ((1 / T) - (1 / stdTemp))) / conversionC
        HenryConstListHONO = T * KHHONO * np.exp(enthalpyHONO * ((1 / T) - (1 / stdTemp))) / conversionC
        HenryConstListNO = T * KHNO * np.exp(enthalpyNO * ((1 / T) - (1 / stdTemp))) / conversionC
        HenryConstListN2O = T * KHN2O * np.exp(enthalpyN2O * ((1 / T) - (1 / stdTemp))) / conversionC

        HenryConstList = dict(O2=HenryConstListO2, CO2=HenryConstListCO2, NH3=HenryConstListNH3, HONO=HenryConstListHONO, NO=HenryConstListNO,N2O=HenryConstListN2O)
    else:
        HenryConstList = dict(O2=HenryConstListO2, CO2=HenryConstListCO2)

    Density_air = air_density(Tdist,0,pa) * 1000

    if len(np.shape(Tdist)) == 0:
        pKList = np.zeros((1,1,9))
    else:
        pKList = np.zeros((np.shape(Tdist)[0],np.shape(Tdist)[1],9))

    pKList[:,:,0] = 2835.76 / T - 0.6322 + 0.001225 * T #ammonia dissociation
    pKList[:,:,1] = 3404.71 / T - 14.8435 + 0.032786 * T #HCO3- CO2 dissociation
    pKList[:,:,2] = 2902.39 / T - 6.498 + 0.02379 * T #HCO3- CO3 2- dissociation
    pKList[:,:,3] = -1228.732 - 0.299444 * T + 35512.75 / T + 485.818 * np.log10(T) #calciate complexation
    pKList[:,:,4] = 171.9065 + 0.077993 * T - 2839.319 / T - 71.595 * np.log10(T) # calciate precipitation 
    pKList[:,:,5] = - np.log10(KaHONO * np.exp(enthalpyHONOoverR * ((1 / T) - (1 / stdTemp)))) # HONO dissociation
    pKList[:,:,6] = -68.2401 +3221.51 / T +25.0627 * np.log10(T)  # gypsum dissolution
    pKList[:,:,7] = 1.24 +0.0036* T # CaSO40 complexation
    pKList[:,:,8] = -1*np.log10(np.exp(140.932 - 13446/T - 22.4773*np.log(T))) # water activity
    #https://aiche.onlinelibrary.wiley.com/doi/pdfdirect/10.1002/aic.690240605?casa_token=FHy3tdvNfR4AAAAA:Dl1s0JZ0mOL-MWcooe12GkKnJCJL6OjEsgHkeEbniRETU1ytNH7aIcFstgvb-DNaTcztwXRuglHtyK8
    return MumaxT,HenryConstList,pKList,Density_air


def chem_gypsum_calcite_aw_Z_reaction(t,y, Z, ions = None, pKs = None, hlist = None, gclist = None, rlist= None, kla=None, gasflow = None, theta = None,totalw=None, SA_gypsum=None, SA_calcite=None,fCA=1, gypsum_mod=True,calcite_mod=True): 
    
    # y0 = 'O2'
    # y1 = 'CO2'
    # y2 = 'HCO3-'
    # y3 = 'CH2O'
    # y4 = 'NO3-'
    # y5 = 'NH3'
    # y6 = 'NO2-'
    # y7 = 'HONO'
    # y8 = 'NO'
    # y9 = 'N2O'
    # y10 = 'NH4+'
    # y11 = 'CO3 2-'
    # y12 = 'Ca2+'
    # y13 = 'CaCO3*'
    # y14 = 'CaCO3'
    # y15 = 'H+'
    # y16 = 'CaSO42H2O'
    # y17 = 'SO42-'
    # y18 = 'CaSO4*'

    # y19 = 'O2' (g)
    # y20 = 'CO2' (g)
    # y21 = 'NH3' (g)
    # y22 = 'HONO'(g) 
    # y23 = 'NO' (g) 
    # y24 = 'N2O' (g)

    # y25 = 'O2' efflux (g) 
    # y26 = 'CO2' efflux (g)
    # y27 = 'NH3' efflux (g)
    # y28 = 'HONO' efflux (g) 
    # y39 = 'NO' efflux (g) 
    # y30 = 'N2O' efflux (g)
    
    aw = 1- 0.017*(np.sum(y[:len(ions)])+Z-y[14]-y[16])
    # Algebraic Equations : ionic interaction/tempearature dependency is not included in this case
    Kaw = 10 ** (- pKs[8])
    #y = np.real(y)
    #y[y <= 0] = 1e-16

    totIon = ions.dot(y[:len(ions)]) - y[15] + Z
    y11 = 0.5 * (np.sqrt(totIon * totIon + 4 * Kaw) - totIon)
    
    if y11 < Kaw:
        frac = np.sqrt(totIon * totIon + 4 * Kaw)
    else:
        frac = Kaw / y11
    
    
    I = 0.5 * ((np.multiply(ions,ions)).dot(y[:len(ions)]) +Z + frac)
    absT = 273.15 + theta
    epsilon = 87.74 - 0.40008 * theta + 0.0009398 * theta ** 2 - 1.41e-06 * theta ** 3
    B = 50.3 / np.sqrt(epsilon * absT)
    G = 1825000.0 * (epsilon * absT) ** (- 1.5)
    # # NH4+ HCO3- CO3^2- H+ Ca2+ NO2-
    # a = np.array([2.5,4,4.5,9,6,3])
    # ionicValency = np.array([1,- 1,- 2,1,2,- 1])
    # NH4+, HCO3- CO3^2- H+ Ca2+ SO42- OH- NO2-
    a = np.array([2.5, 5.4,5.4,9,5,5,3.5, 3]) #From Truesdell and Jones (1974)
    ionicValency = np.array([1, -1,-2,1,2,-2,-1, -1])
    
    #extended Debye-Heckel
    pGamma = np.multiply(G,ionicValency ** 2.0) * (np.sqrt(I) / (1 + np.multiply(np.multiply(a,B),np.sqrt(I))))
    #pGamma =zeros(6,1);
    
    Ka = 10 ** (- pKs[0] + pGamma[3] - pGamma[0])
    K1c = 10 ** (- pKs[1] + pGamma[3] + pGamma[1])
    K2c = 10 ** (- pKs[2] + pGamma[3] + pGamma[2] - pGamma[1])
    KcaComp = 10 ** (- pKs[3] + pGamma[2] + pGamma[4])
    KcaPrecip = 10 ** (- pKs[4] + pGamma[2] + pGamma[4])
    Kahono = 10 ** (- pKs[5] + pGamma[7] + pGamma[3])
    
    KgypsumDiss = 10 ** (- pKs[6] + pGamma[4] + pGamma[5])
    KgypsumComp = 10 ** (- pKs[7] + pGamma[4] + pGamma[5])


    
    # #k0 = 10 ** 10 / (24 * 60 * 60)
    # k0 = np.exp(1246.98 - (6.19*10**4)/absT-183.0*np.log(absT))
    # # reaction for H2CO3* - CO2 hydration CO2eq = KH*CO2(g) :: dissolbed CO2 + H2CO3    
    # reac0 = k0*(y[1]-aw*y[21]*hlist[1]) 

    # CO2 + H2O -> H+ + HCO3-     
    # k1 = np.exp(1246.98 - (6.19*10**4)/absT-183.0*np.log(absT))
    # S = 1000*I/(19.920-1.0049*I) # salinity estimation (Millero, F. J. The marine inorganic carbon cycle. Chemical reviews 2007, 107 (2), 308-341)
    # A4 = 499002.24*np.exp(4.2986*10**-4*S**2 + S*5.75499*10**-5)
    # k2 = A4*np.exp(-90166.83/(8.31*absT))/Kaw
    #k1 = 2221 / (24 * 60 * 60)    
    #k2 = 7.19*10**8/(24*60*60) #/s
    #reac1 = (k1 + k2*frac) * (y[1] - y[2] * y[15] /(aw*K1c))
    
    reac1 = (0.02570*fCA + 8321.76*frac) * (y[1] - y[2] * y[15] /(aw*K1c))
    reac2 = 0
    
    # HCO3 - -> H+ + CO3 2-
    #k3 = 10 ** 10 / (24 * 60 * 60)
    #reac3 = k3 * (y[2] - y[15] * y[11] / K2c)
    
    reac3 = 115740.74 * (y[2] - y[15] * y[11] / K2c)

    # NH4+ --> H+ + NH3
    #k4 = 10 ** 10 / (24 * 60 * 60)
    #reac4 = k4 * (y[10] - y[15] * y[5] / Ka)
    reac4 = 115740.74 * (y[10] - y[15] * y[5] / Ka)
    
    # Nitrous acid and nitrite reaction
    #k5 = 10 ** 10 / (24 * 60 * 60)
    #reac5 = k5 * (y[7] - y[15] * y[6] / Kahono)
    reac5 = 115740.74 * (y[7] - y[15] * y[6] / Kahono)
    
    # complexation of calcite
    #kcom = 10 ** 10 / (24 * 60 * 60)
    #reac6 = kcom * (y[13] - KcaComp * y[11] * y[12])
    reac6 = 115740.74 * (y[13] - KcaComp * y[11] * y[12])

    # Precipitation of calcite
    SIcalcite = y[11] * y[12] / KcaPrecip

    if not calcite_mod:
        #kpre = 10 ** 10 / (24 * 60 * 60)
        #reac7 = kpre * (y[14]/totalw - SIcalcite)
        reac7 = 115740.74 * (y[14]/totalw - SIcalcite)
    else:
        #Hs = y11*pGamma[3]            
        if SIcalcite > 1:
            k_preci = SA_calcite*10**-7
            #k_preci = SA_calcite*calcite_precip_rate(absT, aw, y[13], Hs, y11*10 **(-pGamma[3]), y[1], y[12]*10 ** (-pGamma[4]), y[2]*10 **(-pGamma[1]), y[11]*10**(-pGamma[2]), 10 ** (- pKs[4]), 10 ** (- pKs[2]))
            reac7 = -1*k_preci*(SIcalcite-1)
        else:
            if y[14] >= 0:
                k_diss = SA_calcite*10**-5
                #k_diss = SA_calcite*calcite_diss_rate(absT, aw, y[13], Hs, y11*10 **(-pGamma[3]), y[1], y[12]*10 ** (-pGamma[4]), y[2]*10 **(-pGamma[1]), y[11]*10**(-pGamma[2]), 10 ** (- pKs[4]), 10 ** (- pKs[2]))
                reac7 = k_diss*(1-SIcalcite)
            else:
                reac7 = 0
    
        
    # Autoxidation of NO to NO2- 
    
    #kno = 9e6
    #reac8 = kno*y[8]*y[8]*y[0]
    reac8 = 9e6*y[8]*y[8]*y[0]

    #reac8 = 0 
    # 3HNO2 -> 2NO + NO3- + H+ * H2O
    #khno2 = 1.34e-6
    #reac9 = khno2*(y[7]**4)/y[8]
    reac9 = 0
    
    # push pH to balance change neutrality 
    #k10 = 10 ** 10 / (24 * 60 * 60)
    #reac10 = k10 * (y[15] - y11)
    reac10 = 115740.74 * (y[15] - y11)

    
    # complexation of CaSO4*
    #kpre = 10 ** 10 / (24 * 60 * 60)
    #reac11 = kpre * (y[19] - y[12]*y[18]/KgypsumComp)
    reac11 =  115740.74 * (y[18] - y[12]*y[17]/KgypsumComp)


    # dissolution of gypsum
    if gypsum_mod:
        '''
        For instance, Jeschke et al. (2001)
        reported a gypsum dissolution rate constant k of 1.3 × 10-3 mol/m2/s at near to neutral
        pH, whereas Colombani (2008) reported values ranging between 3 and 7 × 10-5
        mol/m2/s based also on a review of previous work. A recent study by Feng et al. (2017)
        on the dissolution of gypsum (010) cleavage surface by digital holographic microscopy
        reported an average dissolution flux of 3.0 × 10-6 mol/m2/s in deionized water
        
        '''    
        # SA_gypsum #model specific areea for gypsum resolution (0.12m^2/g)
        # Temperature dependence with the activation energe 34 kJ/mol
        
        SIgypsum = aw*aw*y[12]*y[17]/KgypsumDiss
        if SIgypsum > 1:
            kpre = SA_gypsum*2.5*10**-6 # Reiss et al (2021) Minerals 2021, 11(2), 141; https://doi.org/10.3390/min11020141 
            reac12 = -1*kpre*(SIgypsum**0.5-1)**2 #precipitates 
        else:
            if y[16] >= 0:
                k_diss = SA_gypsum*1.3*10**-3
                reac12 = k_diss*(1-SIgypsum) #dissolution of gypsum
            else:
                reac12 = 0
        
    else:
        # Precipitayion of gypsum
        #kpre = 10 ** 10 / (24 * 60 * 60)
        #reac12 = kpre * (y[17]/totalw - y[12] * y[18]/KgypsumDiss)
        reac12 = 115740.74 * (y[16]/totalw - y[12] * y[17]/KgypsumDiss)
    
    
    #gas transfer 
    kLO2 = kla*hlist[0]
    o2trans = kLO2*(hlist[0]*y[19]-y[0]) 
 
    kLCO2 = kla*hlist[1]
    co2trans = kLCO2*(hlist[1]*y[20]-y[1]) 
    
    #volatile gases 
    kLNH3 = kla/hlist[2]
    nh3trans = kLNH3*(hlist[2]*y[21]-y[5]) 
    #volatile gases 
    kLHONO = kla/hlist[3]
    honotrans = kLHONO*(hlist[3]*y[22]-y[7]) 
    
    kLNO = kla*hlist[4]
    notrans = kLNO*(hlist[4]*y[23]-y[8]) 

    kLN2O = kla*hlist[5]
    n2otrans = kLN2O*(hlist[5]*y[24]-y[9]) 
    
    # gasflow input 
    
    o2input = gasflow*(y[19]-gclist[0]) 
    co2input = gasflow*(y[20]-gclist[1]) 
    nh3input = gasflow*(y[21]-gclist[2]) 
    honoinput = gasflow*(y[22]-gclist[3]) 
    noinput = gasflow*(y[23]-gclist[4]) 
    n2oinput = gasflow*(y[24]-gclist[5]) 
    
    
    dy0 = -reac8 +o2trans + rlist[0]# O2
    dy1 = - reac1 - reac2 + co2trans + rlist[1] # CO2
    dy2 = reac1 + reac2 - reac3 + rlist[2]#HCO3- 
    dy3 = + rlist[3]#CH2O
    dy4 = reac9/3 + rlist[4]# NO3- 
    dy5 = reac4 + nh3trans + rlist[5]#NH3
    dy6 = reac5 + reac8 + rlist[6]#NO2- 
    dy7 = - reac5 - reac9 +honotrans + rlist[7]# HONO
    dy8 = -reac8 + 2*reac9/3 + notrans + rlist[8]# NO 
    dy9 = +n2otrans + rlist[9]#N2O 
    dy10 = -reac4 + rlist[10]#NH4+
    dy11 = reac3 + reac6 + reac7/totalw + rlist[11]#CO3 2-   
    dy12 = reac6 + (reac7 + reac12)/totalw +reac11 + rlist[12]# Ca2+
    dy13 = - reac6# CaCO3*
    dy14 = - reac7 #CaCO3 in mg
    dy15 = -reac10#H+

    dy16 = -reac12#'CaSO42H2O' in mg
    dy17 = reac11+reac12/totalw#'SO42-'
    dy18 = -reac11#'CaSO4*'
    
    dy19 = -o2trans -o2input #'O2' (g)
    dy20 = -co2trans - co2input#'CO2' (g)
    dy21 = -nh3trans - nh3input#'NH3' (g)
    dy22 = -honotrans - honoinput #'HONO'(g) 
    dy23 = -notrans - noinput#'NO' (g) 
    dy24 = -n2otrans - n2oinput#'N2O' (g)
    
    # efflux addition
    dy25 = -o2trans # 'O2' efflux (g)
    dy26 = -co2trans # 'CO2' efflux (g)
    dy27 = -nh3trans# 'NH3' efflux (g)
    dy28 = -honotrans #  'HONO' efflux(g) 
    dy29 = -notrans #  'NO' efflux (g) 
    dy30 = -n2otrans # 'N2O' efflux (g)

    return [dy0, dy1, dy2, dy3, dy4, dy5, dy6, dy7, dy8, dy9, dy10,dy11, dy12, dy13, dy14, dy15, dy16, dy17, dy18, dy19, dy20, dy21, dy22,dy23,dy24,dy25, dy26, dy27, dy28, dy29, dy30]
 


def plot_results(PATH,fname_target,tempList,effluxList,timeConcDist,timeLine,saveFigures=False,initD=18,pltD=7): 
    '''
    This fuction outputs the figure similar to the Figure 7 in the paper Kim et al 2024,
    Model simulates the entire period of the example input data (e.g. for winter from 01-01-2019 to 01-02-2019)
    Setting initD (inital day for plotting) and pltD (the duration for plotting), the figure can be generated for the selected time period
    for instance, setting initD = 18 and pltD =7, will show plots from 18th Jan to 25(18+7)th Jan. 
    CO2 efflux (Fc)- field observations and model predictions, local pH values, and input temperautre    
    '''    
    fig, axs = plt.subplots(3,1, sharex=True,figsize=(5,5),gridspec_kw={'height_ratios': [2,1,1]},dpi=300)
    plt.rcParams['font.size'] = '12'
        
    # get labels for plotting (changing from elapsed hours to the time of the observation)
    tempList = tempList.set_index(pd.to_datetime(tempList.time))
    timeList = tempList.resample('1D').mean(numeric_only=True).reset_index().time.to_list()
    tlabels = [] #lables for x axis in month-day format
    for dt in timeList:
        tlabels.append(dt.strftime("%b-%d"))
    tempList['intp_time'] = (tempList.index - tempList.index[0]).total_seconds()/3600
    tempList = tempList.reset_index(drop=True).set_index('intp_time')
    timeL = np.arange(0,tempList.index.max()+1,1) 
    tempList = tempList.reindex(timeL).interpolate()# Change hourly data to 10mins data        
    tempList.time = pd.to_datetime(tempList.time)
    
    # Mark nighttime (par<5) 
    axs[0].fill_between(tempList.index, 0, 1, where=tempList.par < 5, color='gray', ec = 'gray',alpha=0.2,linewidth=0.0, transform=axs[0].get_xaxis_transform())
    axs[1].fill_between(tempList.index, 0, 1, where=tempList.par < 5, color='gray', ec = 'gray',alpha=0.2,linewidth=0.0, transform=axs[1].get_xaxis_transform())
    axs[2].fill_between(tempList.index, 0, 1, where=tempList.par < 5, color='gray', ec = 'gray',alpha=0.2,linewidth=0.0, transform=axs[2].get_xaxis_transform())
        
    chem = 'CO2'    
    #convF = 1/(180*180*180) #change to m2 units:: should be changed when consider profile
    #convF = 1/(1.8*1000000)#change to m2 units:: should be changed when consider profile
    convF = 0.0082*0.0082*12#change to m2 units:: should be changed when consider profile
    efs = effluxList[chem]*1000000/44.01/((timeLine[1]-timeLine[0])*3600)/convF
    axs[0].plot(tempList.index,efs,c='k',lw =2, ls='-', label='Model prediction')
    axs[0].plot(tempList.index, tempList.fc, marker='o',markersize=5,c='r',mec='k',mew= 0.5, alpha=0.6,ls='',label='observation')
    axs[0].axhline(y=0,c='k',ls='--',lw=0.5)

    # Plot pH    
    axs[1].plot(tempList.index,timeConcDist['pH'].ravel(),lw =2, ls='-',c='k', label='soil')

    # Plot input variables, air temperature, soil (5cm) temperature
    axs[2].plot(tempList.index,tempList.ta,lw=2,alpha=0.8,label='air')
    axs[2].plot(tempList.index,tempList.ts,lw=2,alpha=0.8,label='soil (5cm)')
    
    #set xticks 
    axs[0].set_xticks(np.arange(0,timeLine[-1]+1,24))
    axs[0].set_xticklabels([])
    axs[2].set_xticks(np.arange(24*initD+12,24*(initD+pltD)+12,24))
    axs[2].set_xticklabels(tlabels[initD:initD+pltD], rotation=30)

    #set yticks
    axs[0].set_yticks(np.arange(-1,1.1,0.5))
    axs[1].set_yticks([7.6, 7.8,8.0])
    
    #set xlim 
    axs[0].set_xlim(24*initD, 24*(initD+pltD)-1)
    axs[1].set_xlim(24*initD, 24*(initD+pltD)-1)
    axs[2].set_xlim(24*initD, 24*(initD+pltD)-1)
    #set ylim 
    axs[0].set_ylim(-1, 1)
    axs[1].set_ylim(7.5,8.1)
    axs[2].set_ylim(-2, 60)
    
    #Set grids
    axs[2].grid()
    axs[1].grid()
    axs[0].grid()        
    
    #Set y lables
    axs[0].set_ylabel('F$_c$ [$\mu$mol m$^{-2}$s$^{-1}$]')
    axs[1].set_ylabel(r'pH [-]')
    axs[2].set_ylabel('T [$^{\circ}$C]')

    #Set legends
    axs[0].legend(ncol=2, loc='lower right')
    axs[2].legend(ncol=3)
    
    #Add the title
    plt.suptitle(fname_target)
    plt.tight_layout()
            
    if saveFigures:
        plt.savefig(os.path.join(PATH, "model_output", fname_target+".png"),dpi=300)

def main_tabernas_biocrusts(PATH, crust_name, realTinput, respif, fname_target, filesave=True, plotFigures=True,figuresave=False):
    '''

    Parameters
    ----------
    PATH : str
        PATH to call functions and save data and figures.
    crust_name : str
        Name of biocrust types. Here we provide examples of 
        physical depositional crusts (PD),
        mature cyanobacterial crusts (MC)
        lichen crusts dominated by XXXX (SD)
        Choose one case among ['PD', 'MC', 'SD']
    realTinput : TYPE
        DESCRIPTION.
    respif : TYPE
        DESCRIPTION.
    fname_target : str
        name of the files to save outcome
    filesave : boolean, optional
        True when the outcome should be saved. The default is True.
    plotFigures : TYPE, optional
        DESCRIPTION. The default is True.
    figuresave : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    Aside from the saved npz file, three return values are assigned for plotting
    timeLine : list
        DESCRIPTION.
    effluxList : dict
        DESCRIPTION.
    timeConcDist : dict
        DESCRIPTION.

    '''
    

    # Call function to load soil properties and field observations
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
    plottt = 60 # 60 minutes: 1 hour
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

    # This project takes the simple upscaling strategies instead of modelling vertical profiles for computational reasons
    nx = 1 # spatial element in horizontal direction
    ny = 1 # spatial element in vertical direction (along the soil depth)
    

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

            gravimetic_swc = df_obs['swc_md'].iloc[examineT]/1.8
        
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

            # When the simulation domain is assumed heterogeneous, nx and ny can be larger than 1 
            for xx in range(nx):
                for yy in range(ny):
                    
                    InvasedIsland =1
                    kla = 1.e-9*patchArea/waterVolume
                    gasflow = InvasedIsland*0.01/60
                    pKs = pKList[xx, yy, :]
                    theta = thetaTC
                    Z = sitesN['Zion'][xx, yy]/totalweight[xx, yy]
                    totalw = totalweight[xx, yy]
                    
                    # Conversion of CO2 reaction rate from mumol/m2.sec to mol/g_solution.sec
                    # (0.00822)*2 is the area conversion factor from 1g cube (1/1.8)**(1/3)=0.822cm length
                    
                    reactionCO2 = respif*1.e-6*(0.00822)**2/totalw

                    rlist = []
                    for chem in list(chemlist):

                        if chem in ['CO2']:
                            rlist.append(reactionCO2)
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
                        effluxList[chem][examineT, 0] = effluxList[chem][examineT,0] + soln.y[i, -1]*list_compounds.loc[chem].molwt*totalw  # unit in g/g_soil.sec 
    
            for chem in chemlist_woH:
                sitesC[chem] = sitesCtemp[chem]
            
            sitesC['H+'] = sitesCtemp['H+']
            sitesC['pH'] = -np.log10(sitesCtemp['H+']*0.001)
              
                
            for chem in list(sitesC)[:-1]:
                sitesN[chem] = sitesC[chem]*waterVolume

        # write values for every duration of plottt
        for chem in list(sitesC)[:-1]:
            timeConcDist[chem][:,:,examineT] = sitesC[chem]
            sitesN[chem] = sitesC[chem]*waterVolume
            
        timeConcDist['pH'][:,:,examineT] = sitesC['pH']

    print('The end')
    
    if filesave: # Save outcomes 
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
        print('file saved!')

    return timeLine,effluxList, timeConcDist

# %%

def main():
    
    description = textwrap.dedent("""\   
        This script is to
        ORCID:https://orcid.org/0000-0002-3942-3743
    """)
        
    parser = argparse.ArgumentParser(description=description, epilog='',
            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('crust_name', type=str, help='PD, MC, SD')
    parser.add_argument('season', type=str, help='winter, summer')
    parser.add_argument('respif', type=float, help='where to cut the training data; qualitiy assurance level of 50 or75')
    args = parser.parse_args()

    crust_name = args.crust_name
    realTinput = args.season
    respif = args.respif
    filesave = True
    saveFigures = True
    plotFigures = True
    
    fname_target = 'sfDBM_'+crust_name+'_'+realTinput+'_with_reac_respif_%.4f'%respif 
    print(fname_target)
    
    PATH = "C:/Users/kimmi/Documents/GitHub/dfDBM_STOTEN"

    # Load pre-saved soil property data
    df_soil_r, df_obs, averageT = load_field_data_soil_properties(PATH, crust_name, realTinput)
    
    # run a main code of soi-gas exchange model (the simplified field DBM)
    timeLine,effluxList, timeConcDist = main_tabernas_biocrusts(PATH, crust_name, realTinput,respif, fname_target, filesave=filesave, plotFigures=plotFigures)
    
    if plotFigures:    
        if respif == 0:
            figure_title = crust_name+' ('+realTinput+')'
        else:
            figure_title = crust_name+' ('+realTinput+') when $R_{CO_2}$ = %.2f [$\mu$mol m$^{-2}$s$^{-1}$]'%respif 
        
        plot_results(PATH, figure_title,df_obs,effluxList,timeConcDist,timeLine,saveFigures=saveFigures) 

if __name__ == '__main__':
    main()
# %%
# PATH = "C:/Users/kimmi/Documents/GitHub/dfDBM_STOTEN"

# filesave = True
# saveFigures = True
# plotFigures = True

# crust_name = 'PD' # ['PD', 'MC', 'SD']
# realTinput = 'winter' #['winter', 'summer']
# respif = 0 #provide in the unit of micromol/m2.sec-1

# fname_target = 'sfDBM_'+crust_name+'_'+realTinput+'_with_reac_respif_%.4f'%respif 
# print(fname_target)

# df_soil_r, df_obs, averageT = load_field_data_soil_properties(PATH, crust_name, realTinput)

# timeLine,effluxList, timeConcDist = main_tabernas_biocrusts(PATH, crust_name, realTinput,respif, fname_target, filesave=filesave, plotFigures=plotFigures)

# if plotFigures:    
#     if respif == 0:
#         figure_title = crust_name+' ('+realTinput+')'
#     else:
#         figure_title = crust_name+' ('+realTinput+') when $R_{CO_2}$ = %.2f [$\mu$mol m$^{-2}$s$^{-1}$]'%respif 
    
#     plot_results(PATH, figure_title,df_obs,effluxList,timeConcDist,timeLine,saveFigures=saveFigures) 
