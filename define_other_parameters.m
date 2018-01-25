function [otherpar] = define_other_parameters
% Define other parameters not defined in xpar or sub_energetic (e.g.
% parameters you might want to vary!)

otherpar    = zeros(3,1);
otherpar(1) = 120e-3;   % 120e-3 K_ei; cytosolic and IM potassium-concentration
otherpar(2) = 20e-3;    % 20e-3 Mg_tot; cytosolic and IM Mg-concentration
otherpar(3) = 20e-3;    % 20e-3 Pi_e; cytosolic and IM phosphate-concentration
% Changing Pi_e changes baseline deltaPsi_m substantially
otherpar(4) = 3.8;      % 3.8 K_ADTP_dyn; Cytosolic ATP production
otherpar(5) = 0.63;     % 0.7 x_ATPK; log(2)/3.5 (value from Huber 2011)
otherpar(6) = 1;        %1 % K_ADTP_cons; Cytosolic ATP consumption

end