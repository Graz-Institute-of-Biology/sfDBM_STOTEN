# sfDBM_STOTEN

   <p align='center'><img src="https://github.com/minsughim/sfDBM_STOTEN/blob/main/Graphical_abstract.png" width="80%"/></p>
   
 - This repository contains the simplified field application of DBM (dfDBM), which provides a minimized script for simulating gaseous C and N emissions from soils and biocrusts under dynamic hydration conditions. Here, we demonstrate an example of using the model to explain the nocturnal CO<sub>2</sub> uptake by dryland soils with limited biological activities ([Kim et al, 2024](http://doi.org/[10.1016/j.scitotenv.2024.171751)).
 
 - Compared to the previous DBM models in MATLAB (please refer to other repositories, [DBM-Biogeoscience](https://github.com/minsughim/DBM--Biogeoscience) and [DBM-for-drying-soils](https://github.com/minsughim/DBM-for-drying-soils)), this sfDBM is written in Python 3. Another distinction is that, in this sfDBM, the individual-based model of microbial communities is simply replaced as a sink and source of substrates (namely, population-based). However, the kinetics to determine gas-liquid partitioning of compounds, local pH determination under the assumption of charge neutrality, are updated by including gypsum.
 
 - The script provided utilises the field measurements of environmental variables such as air, surface, and soil temperature, water content, atmospheric CO<sub>2</sub>, as well as measured soil properties such as porosity, soil pH, CaCO<sub>3</sub> contents, etc.  The soil property data is provided as supplemantry material in [Kim et al, 2024](http://doi.org/[10.1016/j.scitotenv.2024.171751) also in this repository (./input_data/data_replicates*.csv) The field data and detailed measurement campaigns can be found in [Lopez-Canfin et al., 2022]( https://doi.org/10.1016/j.scitotenv.2022.153746).

-  In this repository, we provide a small subset of the data to demonstrate the usage of the model (located in the folder ./input_data/time_series*.csv). The complete dataset can be found in the data repository by [Lopez-Canfin et al., 2024](https://doi.org/10.5281/zenodo.10836173) 
 
 ## System requirements

### Dependencies (tested version)
	- python (3.10.13)
		- numpy (1.26.4)
		- scipy (1.11.4)
		- pandas (1.5.3)
		- matplotlib (3.8.3)

## Installation
The script ./main_sfDBM.py includes all the functions. The main function takes three arguments: the name of the crust (here, field data samples for PD and SD are provided), the season (either winter or summer), and the respiration rate (either positive or negative, in the unit of μmol $m^{-2} s^{-1}$).
To run the model, for example for the case of physical depositional crust (PD) under dry/cool conditions without additional respiration (rate = 0.0), simply execute the main code as follows from the command line after navigating to the same directory:
```
C:\> cd path\to\script\
C:\path\to\script\> python main_sfDBM.py PD winter 0.0
```
# Reference
**Minsu Kim**, Clément Lopez-Canfin, Roberto Lázaro, Enrique P. Sánchez-Cañete, and Bettina Weber (2024) Unravelling the main mechanism responsible for the nocturnal CO2 uptake by dryland soils. *Science of the Total Environment* http://doi.org/[10.1016/j.scitotenv.2024.171751

Clément Lopez-Canfin, Enrique P. Sánchez-Cañete, and Roberto Lázaro, (2024). Hourly time series of soil and atmosphere variables at the experimental site of El Cautivo, Tabernas Desert, Almeria, Spain (February 2018 to December 2019) \[Data set\]. In *Science of the Total Environment* (1.0.0, Vol. 824). *Zenodo*. https://doi.org/10.5281/zenodo.10836173

Clément Lopez-Canfin, Roberto Lázaro, and Enrique P. Sánchez-Cañete, (2022). Water vapor adsorption by dry soils: A potential link between the water and carbon cycles. *Science of the Total Environment* 824, 153746. https://doi.org/10.1016/j.scitotenv.2022.153746
