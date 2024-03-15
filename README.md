# sfDBM_STOTEN

   <p align='center'><img src="https://github.com/minsughim/sfDBM_STOTEN/blob/main/Graphical_abstract.png" width="80%"/></p>
 - This repository of the simplified field application of DBM (dfDBM) provides a script to simulate gaseous C and N emissions from soils and biocrusts under dynamic hydration condition. Here, we demonstarte an example usage of the model that explains the nocturnal CO<sub>2</sub> uptake by dryland soils where biological activities are limited (Kim et al 2024). 
 - Compared to the previous DBM models in MATLAB (please refer to other repositories, DBM-Biogeoscience and DBM-for-drying-soils), this sfDBM is written in Python 3. Another distinction is that, in this sfDBM, the individual based on model of microbial communities is simply replaced as sink and source of substrates (namely, population based). However, the kinetics to determine gas-liquid partitioning of compounds, local pH determiniation under the aussumption of charge neutralitiy are updated by including gypsum.
 - The provided script utilises the environmental variables at field such as temperature (air, surface , soil), water content, atmospheric CO2, and measured soil properties like porosity, soil pH, CaCO<sub>3</sub> contents, etc. The field data and detailed measurement campagins can be found at  Lopez-Canfin et al. 2021.
-  In this repository, we provide small subset of the data to demonstrate the usage of the model (see the folder ./input_data/). The data can be downloaded via … 
 
 # Model simulation
the script ./main_sfDBM.py includes all the functions and the main function. It takes  

# Reference
**Minsu Kim**, Clément Lopez-Canfin, Roberto Lázaro, Enrique P. Sánchez-Cañete, and Bettina Weber.  Unravelling the main mechanism responsible for the nocturnal CO2 uptake by dryland soils

