function [time, J_ATP_cyto, J_ATP_mito, J_ATP_tot,...
    J_ATP_cyto_prop, J_ATP_mito_prop] = exptData_ATP 
% Experimental data from Padova WT cells
% From 'ATP calculations Padova 020517.xls'
% Only 2 data points calculated, before and after Oligo

time = [8 25]; 

J_ATP_cyto = [15.5 39.8]; 
J_ATP_mito = [24.0 0.5]; 
J_ATP_tot = [43.1 40.3]; 

J_ATP_cyto_prop = [39 99];
J_ATP_mito_prop = [61 1]; 

end