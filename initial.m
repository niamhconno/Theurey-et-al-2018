function [Init]=initial(ADP_e, ATP_e, Ctot0, Qtot0)
%
%    Sets the initial conditions for the state variables
%
Init(1)  = 10^(-8);   % Matrix hydrogen concentration pH=8              H_x
Init(2)  = 0.050;     % Matrix potassium concentration [50 mM]          K_x
Init(3)  = 0.001;     % Matrix magnesium concentration [1 mM]           Mg_x
Init(4)  = 397e-6;    % Reduced NADH in matrix [1.5mM]            NADH_x
Init(5)  = 0.0008;    % Reduced ubiquinol in matrix [0.8 mM]            QH2
Init(6)  = 1.0e-3;    % Reduced cyt-c in IMS [1 mM]                     Cred
Init(7)  = 2.6e-5;    % O2 unused as state variable, parameter only
Init(8)  = ATP_e;     % Matrix free ATP (set to external as initial)    ATP_x
Init(9)  = ADP_e;     % Matrix free ADP (set to external as initial)    ADP_x
Init(10) = 778e-6;         % Matrix ATP bound to magnesium              ATP_mx
Init(11) = 200e-6;         % Matrix ADP bound to magnesium              ADP_mx
Init(12) = 10e-3;     % Matrix inorganic phosphate [10 mM]              Pi_x
Init(13) = 4.6e-3;    % IMS free ATP                                    ATP_i
Init(14) = 230e-6;    % IMS free ADP                                    ADP_i
Init(15) = 0;         % IMS free AMP                                    AMP_i
Init(16) = 4.6e-3;    % IMS ATP bound to magnesium
Init(17) = 226e-6;    % IMS ADP bound to magnesium
Init(18) = 10e-3;     % IMS inorganic phosphate [10 mM]                 Pi_i
Init(19) = 150;       % Mitochondrial membrane potential [150 mV]       dPsi
Init(20) = Ctot0;     % total cyt-c in IMS                              Ctot0
Init(21) = Qtot0;     % total ubiquinol                                 Qtot0
Init(22) = 10^(-8.0); % IMS hydrogen concentration pH=8                 H_i
Init(23) = ATP_e;     % cytosolic ATP (given according text)            ATP_c
Init(24) = ADP_e;     % cytosolic ADP (given according text)            ADP_c

end