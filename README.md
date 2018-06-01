# Theurey-et-al-2018 Aging Cell
Code to simulate model and generate data plotted in Figures 1B-D; 2D,F; 4A

Code created in MATLAB R2017a. 
If there are any issues with this code please contact niamhmconnolly@rcsi.ie

To generate data for all Figures listed above
Run Beard_NC_simulatePopulation.m to simulate population of cells
User input required
% Line 72: Set number of simulations (currently set to 10)
% Line 86: Set what experiment to simulate (currently set to 3 = standard Seahorse OCR)
% Line 96: Set what putative transgenic defects to simulate (all currently set to 1 = wildtype)
Additionally, a break-point should be set at end of code to extract data

To explicitly reproduce Fig 1B:
Run PlotPopulationwLitValues.m

To explicitly reproduce Figs 1C and 2F
Run plotExptBoxplot.m (individual sections for each figure)
(Note: 2018-June-01 - data needs to be added to this code for it to be standalone)
