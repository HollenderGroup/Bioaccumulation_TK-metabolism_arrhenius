# Bioaccumulation_TK-metabolism_arrhenius
This model integrates time varrying exposure (chemical) with phase one biotransformation into a one-compartment TK model. 
A vaildation script for fluctuating temperature (i.e. mescosm study)- using the Arrhenius Equation is provided.

This project repository is connected to the journal article:
‘From laboratory to mesocosm: Temperature and species dependent toxicokinetics of the fluorinated herbicide Diflufenican[’ by Raths et al. 2026 (doi TBA).

The key element of this repository is the R-integration of a one-compartment TK model including phase-one biotransformation (TK-metabolism) and dynamic exposure previously developed in Matlab using the BYOM platform (https://www.debtox.info/byom.html).

How to run this model for your own work?
1. Make sure to provide all necessary input data (tissue and water cocnentrations), water temperatures
2. Open "TK_metabolism_arrhenius_calibration.R" and follow the instructions. Do not forget to update the molar weights of your analytes!
3. For validation use "TK_metabolism_validation_mesocosm.R" and follow the instructions there. An input file including the Arrhenius fits (*Arrhenius_results.xlsx") may be needed.

