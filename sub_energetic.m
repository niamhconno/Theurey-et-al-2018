function [f] = sub_energetic(t,x,xpar,otherpar,time,oligo,rotenone,AA,CIV,FCCP,energy)
%==================================================================
% Code supplement to Huber, Duessmann, Rehm and Prehn
%
% Bioenergetic part of the model including cytosolic extensions
% Heinrich Huber, Royal College of Surgeons in Ireland, 2010
%
% Email comments and suggestions to heinhuber@rcsi.ie
%==================================================================
%==================================================================
% The bioenergetic model is derived from the great work of
% Daniel A. Beard, Department of Physiology, Medical College of Wisonsin,
% Milwaukee, WI.
%
% Beard DA. (2005) A Biophysical Model of the Mitochenondrial Respiratory
% System and Oxidative Phosphorylation. PLoS Comp Bio. 1 (4) e36
%
% Including the ATP/ADP transfer modelling from:
% Korzeniewski B. (1998), Biophys Chem 75(1), 73-80
%==========================================================================
% We have included the following features
%     - Facultative blockage of ATP synthase
%     - Cytosolic ATP production and consumption
%
%==========================================================================
%==========================================================================
%
% Notes
%     - For consistency and comparability to Beard's original model
%       we kept the constants in seconds and then rescaled them to minutes
%     - To maintain redox equilibrium, always the total and one redox
%       fraction of cyt-c/ubiqinon is affected from release/cleavage
%     - Complex II is not explicitetly considered (similar to Beard).
%       Caspase-3 has been reported to cleave both, complex I and II.
%       Likewise, as they work in parallel, they are lumped together
%==========================================================================
%
% INPUTS:
%   t - time
%   x - vector of state variables
%
% OUTPUTS:
%   f - vector of time derivatives of concentrations
%
% Global INPUT (*No longer applicable without apoptosis sub-function*)
%   pp_casp3:  Spline parameter of temporal concentrations for free, active
%              caspase-3 (in uM) starting from time of depolarization (i.e.
%              shifted)
% Global OUTPUT
%   C - Output vector for plotted concentration
%
%========================================================================================
%
global C pp_casp3
               
%*************************************************************************
% Concentration state variables: (rename them for clarity)
%*************************************************************************
H_x    = x(1);                      % Matrix hydrogen concentration
K_x    = x(2);                      % Matrix potassium concentration
Mg_x   = x(3);                      % Matrix magnesium concentration
NADH_x = x(4);                      % Reduced NADH in matrix
QH2    = x(5);                      % Reduced ubiquinol in matrix
Cred   = x(6);                      % Reduced cyt-c in IMS
O2     = x(7);                      % O2 unused as state variable, parameter only
ATP_x  = x(8);                      % Matrix free ATP (set to external as initial)
ADP_x  = x(9);                      % Matrix free ADP (set to external as initial)
ATP_mx = x(10);                     % Matrix ATP bound to magnesium
ADP_mx = x(11);                     % Matrix ADP bound to magnesium
Pi_x   = x(12);                     % Matrix inorganic phosphate
ATP_i  = x(13);                     % IMS free ATP
ADP_i  = x(14);                     % IMS free ADP
AMP_i  = x(15);                     % IMS free AMP
ATP_mi = x(16);                     % IMS ATP bound to magnesium
ADP_mi = x(17);                     % IMS ADP bound to magnesium
Pi_i   = x(18);                     % IMS inorganic phosphate
dPsi   = x(19);                     % Mitochondrial membrane potential
Ctot   = x(20);                     % total cyt-c in IMS
Qtot   = x(21);                     % total ubiquinol
H_i    = x(22);                     % IMS hydrogen concentration
ATP_e  = x(23);                     % cytosolic ATP
ADP_e  = x(24);                     % cytosolic ADP
AMP_e  = 0;                         % cytosolic AMP concentration (Molar), Nicholls and others

%*************************************************************************
%   Model parameters (fixed parameters)
%*************************************************************************
%   Thermodynamic parameters
F       = 0.096484;                  % kJ mol^{-1} mV^{-1}
R       = 8314e-6;                   % Universal gas constant (kJ * mol^{-1} * K^{-1})
T       = (273.15 + 37);             % Temperature (K), 37 degree
RT      = R*T;                       % kJ  mol^{-1}

%   Parameters for literature sources, please refer to Beard 2005 PLoS Comp Bio. 1 (4) e36)
dG_C1o  = -69.37;                    % kJ mol^{-1}
dG_C3o  = -32.53;                    % kJ mol^{-1}
dG_C4o  = -122.94;                   % kJ mol^{-1}
dG_F1o  = 36.03;                     % kJ mol^{-1}
n_A     = 3.0;                       % numbers of proteins used by ATP synthase

%   Cytosolic Ion concentration and pH
pH_e    = 7.4;                       % External pH (cytosol)
H_e     = 10^(-pH_e);                % cytosolic hydrogen concentration (Molar)
K_ei    = otherpar(1); %120.0e-3;    % cytosolic and IM potassium-concentration (Nicholls, Bioenergetics 2)
Mg_tot  = otherpar(2); %0.020;       % cytosolic and IM magnesium-concentration (Nicholls, Bioenergetics 2).
Pi_e    = otherpar(3); %0.02;        % cytosolic and IM phosphate-concentration
                           
%if t> time.ATPcyto                    
% This does nothing as commands outside for loop (near beginning 
% of code) are the same as for else commands
%else                                % If fixing ATP (and for ss calculations)
%    ADTP_tot= xpar(4);              % total Adenosine phosphates in cytosol, calculated from Evans_77
%    ATP_e   = state_fact*ADTP_tot;  % State 3 ATP level
%    ADP_e   = ADTP_tot-ATP_e;       % State 3 ADP level
%end

c_inc   = 1;                     % Cyt-c after release set to 0.1%
V_mito  = 0.06;                  % 6% mitochondrial fraction Ward 2007
W_c     = 1/V_mito;              % Volume fraction cytosol/mitochondria
        
NADtot      = xpar(1);           % Calculated from Evans 77.
Ctot0       = xpar(2);           % Initial Cyt-c in IMS, Cred+Cox (M)
Qtot0       = xpar(3);           % Initial IMS Ubiquinol, Q+QH2 (M)

% Parameters for input function
k_Pi1  =  xpar(5);                 % Dehydrogenase flux input
k_Pi2  =  xpar(6);                 % Dehydrogenase flux input
x_DH   =  xpar(7);                 % Dehydrogenase activity
r_DH   =  xpar(8);                 % Input-flux: Initial disturbance of equilibrium

% Parameters for complex I
x_C1   =  xpar(9);                 % Complex I activity

% Parameters for complex III
x_C3   =  xpar(10);                % Complex III activity
k_Pi3  =  xpar(11);                % Complex III / Pi parameter
k_Pi4  =  xpar(12);                % Complex III / Pi parameter

% Parameters for complex IV
x_C4   =  xpar(13);                % Complex IV activity
k_O2    = xpar(14);                % kinetic constant for complex IV (M)

% Parameters for ATP synthase and Mg-binding to ATP?ADP
x_F1   =  xpar(15);                % F1F0 ATPase activity
K_DT    = xpar(16);                % Mg/ATP binding constant M
K_DD    = xpar(17);                % Mg/ADP binding constant M
x_MgA  =  xpar(18);                % Mg2+ binding activity

% Parameters for Adenosine Transferase
x_ANT  =  xpar(19);                % ANT activity
   % Simulating increased ATP export from the mitochondria)
   % Initially used to model increased energy demand
   %if t > energy.t, x_ANT = x_ANT * 10; end
k_mADP =  xpar(20);                % ANT Michaelis-Menten constant

% Proton leaks
x_Hle  =  xpar(21);                % Proton leak activity

% Parameters for OM transporters
x_Ht    = xpar(22);                % mito outer membrane permeability to protons (micron s^{-1})
gamma   = xpar(23);                % mito outer membrane area per unit volume micron^{-1}
x_A     = xpar(24);                % mito outer membrane permeability to nucleotides
x_Pi2   = xpar(25);                % mito outer membrane permeability to phosphate

% Phosphate-Hydrogen-Cotransport
k_dHPi  = xpar(26);                % H/Pi binding constant Molar
k_PiH   = xpar(27);                % H+ / Pi ? co-transport Michaelis constant
x_Pi1   = xpar(28);                % H+ / Pi ? co-transport activity

% Potassium-Hydrogen-Antiport
x_KH   =  xpar(29);                % K+ / H+ antiporter activity

% Matrix buffer and membrane capacitance
x_buff =  xpar(30);                % Inner Matrix hydrogen buffer capacity
CIM     = xpar(31);                % Inner Membrane capacitance

% Release and compartment parameters
W_m     = 0.143/0.20;               % mitochondrial water space (ml water / ml mito)
W_x     = 0.9*W_m;                  % Matrix water space (ml water / ml mito)
W_i     = 0.1*W_m;                  % IM water space (ml water / ml mito)
t12cyto = 90;                       % cyt-c release half time Waterhouse 2000 (seconds)
Cfin    = Ctot0*W_i/W_c*c_inc;      % final cytochrome c after release in IMS given by re-equilibration of the
                                    % IMS water space with the rest of the cell (0.1%)
                                    % c_inc is used to modify OXPHOS-contributing cyt-c after release

        %*************************************************************************
        % potassium uniporter and adenylate kinase neglected
        %*************************************************************************
        x_K    = 0;                         % Passive potassium transporter activity
        x_AK   = 0e6;                       % AK activity
        K_AK    = 0;                        % Adenelyte Kinase switched off

        %*************************************************************************
        % Balancing moeities
        %*************************************************************************
        %
        NAD_x  = NADtot - NADH_x;
        Q      = Qtot   - QH2;
        Cox    = Ctot   - Cred;
        null = 0;
        %
        %*************************************************************************
        % Other concentrations computed from the state variables
        %*************************************************************************
        %
        ATP_fx  = ATP_x - ATP_mx;
        ADP_fx  = ADP_x - ADP_mx;
        ATP_fi  = ATP_i - ATP_mi;
        ADP_fi  = ADP_i - ADP_mi;
        %
        %*************************************************************************
        % Preventing numerical instabilities
        % Sets to 0 if value < 0
        %*************************************************************************
        %
        Q       = Q*(Q>0);
        QH2     = QH2*(QH2>0);
        Qtot    = Qtot*(Qtot>0);
        Cox     = Cox*(Cox>0);
        Cred    = Cred*(Cred>0);
        Ctot    = Ctot*(Ctot>0);
        ATP_fx  = ATP_fx*(ATP_fx>0);
        ATP_fi  = ATP_fi*(ATP_fi>0);
        ADP_fx 	= ADP_fx*(ADP_fx>0);
        ADP_fi  = ADP_fi*(ADP_fi>0);
        %
        %*************************************************************************
        % ADP/Mg/K binding in E space
        %*************************************************************************
        ADP_me  = ( (K_DD+ADP_e+Mg_tot) - sqrt((K_DD+ADP_e+Mg_tot)^2-4*(Mg_tot*ADP_e)) )/2;
        ADP_fe  = ADP_e - ADP_me;
        Mg_e    = Mg_tot - ADP_me;
        Mg_i    = Mg_e;                         % Mg,K in IM space same as external/cytosolic
        K_i     = K_ei; 
        %
        %*************************************************************************
        % Modelling ATP synthase block (Oligomycin)
        % Previous attempts - see further down for actual ATP synthase
        % block
        %*************************************************************************
        %if t>oligo.t
        %    x_F1 = 0.4*x_F1;      % reduced ATP synthase activity
        %end
        
        %*************************************************************************
        % Calculating Membrane proton motive force and respiration fluxes
        %*************************************************************************
        dG_H    = F*dPsi + 1*RT*log(H_i/H_x);   % Protomotive force
        dG_C1op = dG_C1o - 1*RT*log(H_x/1e-7);
        dG_C3op = dG_C3o + 2*RT*log(H_x/1e-7);
        dG_C4op = dG_C4o - 2*RT*log(H_x/1e-7);
        dG_F1op = dG_F1o - 1*RT*log(H_x/1e-7);
        %
        J_DH    = x_DH*(r_DH*NAD_x-NADH_x)*((1+Pi_x/k_Pi1)/(1+Pi_x/k_Pi2));
        J_C1    = x_C1*(exp(-(dG_C1op+4*dG_H)/RT)*NADH_x*Q - NAD_x*QH2);
        J_C3    = x_C3*((1+Pi_x/k_Pi3)/(1+Pi_x/k_Pi4))*...
            (exp(-(dG_C3op+4*dG_H-2*F*dPsi)/(2*RT))*Cox*QH2^0.5 - Cred*Q^0.5);
        J_C4    = x_C4*(O2/(O2+k_O2))*(Cred/Ctot)*...
            (exp(-(dG_C4op+2*dG_H)/(2*RT))*Cred*(O2^0.25) - Cox*exp(F*dPsi/RT));
        J_F1    = x_F1*(exp(-(dG_F1op-n_A*dG_H)/RT)*(K_DD/K_DT)*ADP_mx*Pi_x - ATP_mx);
        
        %
        %*************************************************************************
        % Modelling ATP transferase Korzeniewski 1998
        %*************************************************************************
        Psi_x = -0.65*dPsi;
        Psi_i = +0.35*dPsi;
        if (ADP_fi > null) || (ATP_fi > null)
            J_ANT  = x_ANT*( ADP_fi/(ADP_fi+ATP_fi*exp(-F*Psi_i/RT))...
                - ADP_fx/(ADP_fx+ATP_fx*exp(-F*Psi_x/RT)) )*(ADP_fi/(ADP_fi+k_mADP));
        else
            J_ANT  = 0;
        end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Modelling facultative caspase-3 feedback to Complex I/II, from Huber 2011
%%% (Caspase-3 activation cleaves Complex I/II)
%%% Set k_cleave 0 to remove this feedback (NB. Apoptosis sub-routine has been
%%% removed so this will not work unless another executable MATLAB file
%%% containing the apoptosis sub-routine has been executed (to output
%%% ppcasp3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k_cleave        = 0;         % cleavage constant       
%k_cleave        = 12/60*1.e6;       % Caspase 3 activity from Stennicke 2000 (in M per seconds)

if (t >= time.release)
  if k_cleave>0
    J_CASP_Qtot  = k_cleave*ppval(pp_casp3,t-time.release)*1e-6 * (Qtot-Qtot0/100);
    % assumption: concentration of caspase 3 in IM space equals cytosolic
    J_CASP_QH2   = k_cleave*ppval(pp_casp3,t-time.release)*1e-6 * QH2;
    % assumption: concentration of caspase 3 in IM space
    % normalization of caspase 3 concentration in molar
  else
    J_CASP_Qtot  = 0;   % No cleavage of total ubiquinone
    J_CASP_QH2   = 0;   % No cleavage of reduced ubiquinone
  end
else
  J_CASP_Qtot  = 0;       % No cleavage of total ubiquinone
  J_CASP_QH2   = 0;       % No cleavage of reduced ubiquinone
end
        
        %*************************************************************************
        % Calculating ionic fluxes
        %*************************************************************************
        %
        H2PIi       = Pi_i*H_i/(H_i+k_dHPi); 
        H2PIx       = Pi_x*H_x/(H_x+k_dHPi);
        J_Pi1       = x_Pi1*(H_x*H2PIi - H_i*H2PIx)/(H2PIi+k_PiH);
        J_Hle       = x_Hle*dPsi*(H_i*exp(F*dPsi/RT)-H_x)/(exp(F*dPsi/RT)-1);
        J_KH        = x_KH*( K_i*H_x - K_x*H_i );
        J_K         = x_K*dPsi*(K_i*exp(F*dPsi/RT)-K_x)/(exp(F*dPsi/RT)-1);
        J_AKi       = x_AK*( K_AK*ADP_i*ADP_i - AMP_i*ATP_i );
        J_AMP       = gamma*x_A*(AMP_e-AMP_i);
        J_ADP       = gamma*x_A*(ADP_e-ADP_i);
        J_ATP       = gamma*x_A*(ATP_e-ATP_i);
        J_Pi2       = gamma*x_Pi2*(Pi_e-Pi_i);
        J_Ht        = gamma*x_Ht*(H_e-H_i);
        %
        J_MgATPx    = x_MgA*(ATP_fx*Mg_x-K_DT*ATP_mx);
        J_MgADPx    = x_MgA*(ADP_fx*Mg_x-K_DD*ADP_mx);
        J_MgATPi    = x_MgA*(ATP_fi*Mg_i-K_DT*ATP_mi);
        J_MgADPi    = x_MgA*(ADP_fi*Mg_i-K_DD*ADP_mi);
        %
        %*************************************************************************
        %   Cytosolic Energy balance for single cell model (see Huber 2011)
        %*************************************************************************
        K_ADTP_dyn = otherpar(4);
        % Increase cytosolic ATP consumption to simulate increased ATP
        % demand when t > energy.t
        if t > energy.t
            K_ADTP_cons = energy.factor*otherpar(6);
        else, K_ADTP_cons = otherpar(6);
        end
        
        if  t > time.ATPcyto           % Cytosolic ATP production and consumption
            % x_ATPK varied to establish desired ss equilibrium between
            % ATP:ADP (See x_ATPK and K_ADTP_dyn calcs)
            x_ATPK   = otherpar(5); %log(2)/3.5 (value from Huber 2011)   
        else
            x_ATPK   = 0;
        end
                    
        J_ATPK   = x_ATPK * (K_ADTP_cons*ATP_e- K_ADTP_dyn*ADP_e);
        %
        %*************************************************************************
        % Modelling ATP synthase block (Oligomycin)
        % This previously fixed the flux through CV (i.e. independent of 
        % substrate/upstream respiratory chain activity), rather than inhibiting
        % its activity. But reducing x_F1 (activity factor), results in
        % only transient decrease in flux through CV.
        %*************************************************************************
        if t>oligo.t
            %if oligo.check == 1
            %    oligo.J_F1_pre = J_F1;
            %    oligo.check = 0;
            %end
            J_F1 = oligo.percent * J_F1; % oligo.percent*oligo.J_F1_pre;       
            %J_ANT  = oligo.percent * J_ANT;      % as consequence no ATP/ADP transferase
        end
        
        %*************************************************************************
        % Modelling Complex I block (Rotenone)
        % Modelling Complex III block (Antimycin A)
        % Modelling Complex IV block (e.g. CN)
        % Modelling increased proton leak (FCCP)
        %*************************************************************************
        if t > rotenone.t,  J_C1 = rotenone.percent * J_C1;  end
        if t > AA.t,  J_C3 = AA.percent * J_C3;  end
        if t > CIV.t,  J_C4 = CIV.percent * J_C4;  end
        if t > FCCP.t,  J_Hle = FCCP.factor * J_Hle;  end
        
        %*************************************************************************
        %   Calculating derivatives for next step
        %*************************************************************************
        %
        f(1)  = x_buff*H_x*( +1*J_DH - (4+1)*J_C1 - (4-2)*J_C3 - (2+2)*J_C4...
            + (n_A-1)*J_F1 + 2*J_Pi1 + J_Hle - J_KH )/W_x; % H_x
        f(2)  = (J_KH + J_K)/W_x; % K_x
        f(3)  = (-J_MgATPx - J_MgADPx)/W_x; % Mg_x
        f(4)  = (+J_DH - J_C1)/W_x; % NADH
        % f(5) includes caspase 3 feedback
        % set to (+J_C1 - J_C3)/W_x without this
        f(5)  = (+J_C1 - J_C3)/W_x - J_CASP_QH2*(t>time.release); % QH2
        % f(6) includes loss of reduced cyt-c with half time t12cyto; 
        % Set to (+2*J_C3 - 2*J_C4)/W_i without this
        f(6)  = (+2*J_C3 - 2*J_C4)/W_i- log(2)/t12cyto* (Cred - Cfin)*(t>time.release);  
        f(7)  = 0; % O2 unused
        f(8)  = (+J_F1 - J_ANT)/W_x; % ATP_x
        f(9)  = (-J_F1 + J_ANT)/W_x; % ADP_x
        f(10) = (J_MgATPx)/W_x; % ATP_mx
        f(11) = (J_MgADPx)/W_x; % ADP_mx
        f(12) = (-J_F1 + J_Pi1 )/W_x;  % Pi_x MOD
        f(13) = (+J_ATP + J_ANT + J_AKi )/W_i; % ATP_i
        f(14) = (+J_ADP - J_ANT - 2*J_AKi )/W_i; % ADP_i
        f(15) = (+J_AMP + J_AKi)/W_i; % AMP_i
        f(16) = (J_MgATPi)/W_i; % ATP_mi
        f(17) = (J_MgADPi)/W_i; % ADP_mi
        f(18) = (-J_Pi1 + J_Pi2 )/W_i; % Pi_i
        f(19) = ( 4*J_C1 + 2*J_C3 + 4*J_C4 - n_A*J_F1 - J_ANT - J_Hle  - J_K )/CIM; %MOD
        % f(20), f(21) include cytochrome-c releas4e and Caspase cleavage 
        % Set to 0 without these (which I don't use)
        f(20) = - log(2)/t12cyto*(Ctot - Cfin)*(t>time.release);         % Cytochrome-c loss with half time t12cyto
        f(21) = - J_CASP_Qtot*(t>time.release);                          % Cleavage of complex I/II total
        f(22) = (- J_DH +(4+1)*J_C1 + (4-2)*J_C3 + (2+2)*J_C4 - (n_A-1)*J_F1 ...
            - 2*J_Pi1 - J_Hle + J_KH + J_Ht)/W_i;                   % Change Intermembrane Hydrogen

        if  t>time.ATPvary
            f(23) = (-J_ATP - J_ATPK)/W_c;      % ATP_c; ATP levels established by mitochondria & glycolysis
            f(24) = (-J_ADP + J_ATPK)/W_c;      % ADP_c
        else
            f(23) = 0;                          % Fixed external/cytosolic ATP for e.g. isolated mitochondria
            f(24) = 0;
        end
        %
        %*************************************************************************
        % Scale Time to minutes (parameters in seconds, but output in minutes)
        %*************************************************************************
        %
        f = 60*f';                          % Scale to minutes
        %
%*************************************************************************
% Output variables
%*************************************************************************
C(1)    = J_DH;           % Dehydrogenase flux (input function)
C(2)    = J_C1;           % Respiration/electron flux through Complex I/II
C(3)    = J_C3;           % Respiration/electron flux through Complex III 
C(4)    = J_C4;           % Respiration/electron flux through Complex IV
C(5)    = J_F1;           % ATP (FoF1) synthase flux
C(6)    = J_ANT;          % ATP transfer matrix IMS (ATP transferase inner membrane)
C(7)    = J_ATP;          % ATP transfer IMS cytosol
C(8)    = J_ADP;          % ADP transfer IMS cytosol
C(9)    = H_x;            % Matrix Hydrogen
C(10)   = Ctot;           % Total available cyt-c (IMS)
C(11)   = Qtot;           % Total available Complex1 (Total ubiquinol)

if k_cleave>0
  C(12) = (k_cleave>0)*ppval(pp_casp3,t-time.release)*(t>time.release);           % free active caspase3 in nM
else C(12) = 0;
end

C(13) = dPsi;                   % Mitochondrial Membrane potential
C(14) = dG_H/F;                 % Membrame Protomotive Force
C(15) = J_Hle;                  % Proton leak flux (Hydrogen leaks)
C(16) = ATP_e/(ADP_e+ATP_e);    % Cytosolic (external Medium for isolated mitochondria) ATP ratio
C(17) = -log10(H_x) - pH_e;     % pH difference between matrix and cytosolic pH
        
toth    = (4+1)*J_C1 +(4-2)*J_C3 +(2+2)*J_C4 + J_KH ;           % total proton extrusion
C(18)   = 100*((4+1)*J_C1 + (4-2)*J_C3 + (2+2)*J_C4)/toth;      % proton extrusion by respiration
C(19)   = 100*(J_KH)/toth;              % proton generation by potassium hydrogen
C(20)   = 100*(J_DH)/toth;              % proton consumption by inpt flux
C(21)   = 100*(n_A-1)*J_F1/toth;        % proton consumption by ATP synthase
C(22)   = 100*2*J_Pi1/toth;             % proton consumption by phosphate-hydrogen
C(23)   = 100*J_Hle/toth;               % proton consumption by proton leaks
C(24)   = ATP_e;                        % cytosolic ATP
C(25)   = NAD_x;                        % Matrix NAD
C(26)   = -log10(H_x);                  % Matrix pH
if t > time.ATPcyto% so long as cyto ATP processes are switched on
    C(27)   = x_ATPK * K_ADTP_dyn * ADP_e;  % Cytosolic ATP production
else C(27) = 0;
end
if t > time.ATPvary     % If cytosolic ATP is 'allowed' to vary
    C(28) = x_ATPK * K_ADTP_cons * ATP_e;   % Cytosolic ATP consumption
else C(28) = 0;
end

end