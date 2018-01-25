function [xpar] = define_model_parameters
%*************************************************************************
% Define model parameters
% Updated from Huber
%*************************************************************************
%
% Metabolites
%
xpar    = zeros(31,1);
xpar(1) = 726e-6;     % NADtot    % NADH calculated from Evans_77
xpar(2) = 2.70e-3;                          % Cyt-c after Beard (total IMS Cyt-C, Cred+Cox, molar)
xpar(3) = 1.35e-3;                          % total IMS Ubiquinol, Q+QH2, molar
xpar(4) = 2.6e-3;         %ADTPtot; total Adenosine phosphates in cytosol
                                            
%*************************************************************************
% FIG. 2 FLUX 1: Parameter for input function
%*************************************************************************
xpar(5) = 0.13413e-3;          % k_Pi1      % Dehydrogenase flux input
xpar(6) = 0.677e-3;          % k_Pi2      % Dehydrogenase flux input
xpar(7) = 0.0737;    %0.0918;      % x_DH       % Dehydrogenase activity
xpar(8) =  4.3;    %45.807;      % r_DH       % Input-flux: Initial disturbance of equilibrium
%*************************************************************************
% FIG. 2 FLUX 2: Parameter for complex I
%*************************************************************************
xpar(9) = 1.0200e+003;         % x_C1       % Complex I activity
%*************************************************************************
% FIG. 2 FLUX 3: Parameter for complex III
%*************************************************************************
xpar(10) = 0.2241;              % x_C3      % Complex III activity
xpar(11) = 0.192e-3;          % k_Pi3     % Complex III / Pi parameter
xpar(12) = 25.31e-3;            % k_Pi4     % Complex III / Pi parameter
%*************************************************************************
% FIG. 2 FLUX 4: Parameter for complex IV
%*************************************************************************
xpar(13) = 3.2131e-004;         % x_C4      % Complex IV activity
xpar(14) = 1.2e-4;              % k_O2      % kinetic constant for complex IV (M)
%*************************************************************************
% FIG. 2 FLUX 5: Parameter for ATP synthase and Mg-binding to ATP/ADP
%*************************************************************************
xpar(15) = 6.8294e+003;         % x_F1      % F1Fo ATPase activity
xpar(16) = 24e-6;               % K_DT      % Mg/ATP binding constant (M)
xpar(17) = 347e-6;              % K_DD      % Mg/ADP binding constant (M)
xpar(18) = 1e6;                 % x_MgA     % Mg2+ binding activity
%*************************************************************************
% FIG. 2 FLUX 6: Parameter for Adenosine Transferase
%*************************************************************************
xpar(19) = 0.0020;              % x_ANT     % ANT activity
xpar(20) = 3.50e-6;             % k_mADP    % ANT Michaelis-Menten constant
%*************************************************************************
% FIG. 2 FLUX 7: Proton leak activity
%*************************************************************************
xpar(21) = 150.0;               % x_Hle     % Proton leak activity  
%*************************************************************************
% FIG. 2 FLUX 8: Parameter for OM transporters
%*************************************************************************
xpar(22) = 2e1;                 % x_Ht      % MOM permeability to protons (micron s^{-1})
xpar(23) = 5.99;                % gamma     % MOM area per unit volume micron^{-1}
xpar(24) = 85.0;                % x_A       % MOM permeability to nucleotides (micron sec^{-1}) 
xpar(25) = 327;                 % x_Pi2     % MOM permeability to phosphate (micron sec^{-1})
%*************************************************************************
% FIG. 2 FLUX 9: Phosphate-Hydrogen-Cotransport
%*************************************************************************
xpar(26) = 10^(-6.75);          % k_dHPi    % H/Pi co-transpor Molar (form factor 1) (binding)
xpar(27) = 0.45082e-3;          % k_PiH     % H+/Pi co-transport activity (form factor 2) (Michaelis constant)
xpar(28) = 3.85e5;            % x_Pi1     % H+/Pi? co-transport activity
%*************************************************************************
% FIG. 2 FLUX 10: Potassium-Hydrogen-Antiport
%*************************************************************************
xpar(29) = 2.9802e7;            % x_KH      % K+ / H+ antiporter activity

xpar(30) = 100;                 % x_buff    % Inner Matrix hydrogen buffer capacity
xpar(31) = 6.7568e-6;           % CIM       % Inner Membrane capacitance

end