# Theurey-et-al-2018 Aging Cell
Code to simulate model and generate data plotted in Figures 1B-D; 2Di,F; 4A

If you make use of this code please cite: 

Pierre Theurey, Niamh M. C. Connolly, Ilaria Fortunati, Emy Basso, Susette Lauwen, Camilla Ferrante, Catarina Moreira Pinho, Alvin Joselin, Anna Gioran, Daniele Bano, David S. Park, Maria Ankarcrona, Paola Pizzo, Jochen H. M. Prehn (2019). Systems biology identifies preserved integrity but impaired metabolism of mitochondria due to a glycolytic defect in Alzheimer’s disease neurons. Aging Cell. 2019; e12924. https://doi.org/10.1111/acel.12924

Code created in MATLAB R2017a. 

If there are any issues with this code please contact niamhmconnolly@rcsi.ie


-To generate data for all Figures listed above

--Run Beard_NC_simulatePopulation.m to simulate population of cells

---User input required

---% Line 72: Set number of simulations (currently set to 10)

---% Line 86: Set what experiment to simulate (currently set to 3 = standard Seahorse OCR)

---% Line 96: Set what putative transgenic defects to simulate (all currently set to 1 = wildtype)

--Additionally, a break-point should be set at end of code to extract data

-- Fig 2Di curves can be seen in the Figure 4 automatically generated by code (top graph in 2nd column) (100 sims)


-To explicitly reproduce Fig 1B:

--Run PlotPopulationwLitValues.m


-To explicitly reproduce Figs 1C and 2F

--Copy relevant data from "Pre-generated data (Fig 1c,2F).xlsx" into a variable in MATLAB (called x)

--Run relevant code sections of "plotExptBoxplot.m" (individual sections for each figure)

