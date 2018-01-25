function [time, ROF_glc_WT_FC, AOF_glc_WT_FC, ROF_glc_WT_mV,...
    AOF_glc_WT_mV] = exptData_TMRM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TMRM experimental data from Padova WT
%%% Mean value from all cells
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Time in mins
time = [0 5 10 15 20 25 30 35 40 45 50 55 60 65 70 75];   

% TMRM foldchange data in Glucose conditions, WT cells; response to
% Rotenone and Oligo
ROF_glc_WT_FC = [0.98 1.01 1.04 1 0.87 0.73 0.68 0.63 0.60 0.53 0.50...
    0.47 0.46 0.42 0.39 0.07];

% TMRM foldchange data in Glucose conditions, WT cells; response to
% Antimycin A and Oligo
AOF_glc_WT_FC = [0.89 0.95 1 0.96 0.81 0.72 0.66 0.62 0.60 0.46 0.35...
    0.29 0.25 0.22 0.20 0.07]; 

% TMRM mV data in Glucose conditions, WT cells; response to
% Rotenone and Oligo. Quantification details in TMRMquantROFglc.mat and
% Evernote
ROF_glc_WT_mV = [149.4 150.3 151.1 149.8 145.5 138.7 135.5 132.5...
    130.1 126.1 123.2 120.9 120.1 117.6 116.4 78.9];

% TMRM foldchange data in Glucose conditions, WT cells; response to
% Antimycin A and Oligo Quantification details in TMRMquantAOFglc.mat and
% Evernote
AOF_glc_WT_mV = [148.7 150.5 151.8 150.5 145.5 141.8 138.9 136.6...
    134.9 127.7 119.9 114.7 111.2 107.9 104.3 76.7];

end