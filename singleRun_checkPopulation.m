function Beard_NC_singleRun_checkPopulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code adapted from crosst.m of Huber 2011 (MSB) and 
%%% regulation.m of Huber 2012 (MBS). 
%%% This function simulates population for all experiments
%%% Plots important steady-state parameter values (and literature values)
%%% and timecourses

%%% This doesn't currently plot the timeseries properly (most of them =0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global C  

%%% Define model parameters
xpar = define_model_parameters;
otherpar = define_other_parameters;

%%% Define simulation time-frames and other time settings
printRequest = 0;   %Set to 1 if you want to print info on time settings
[t_prior, t_start,t_final,t_no_time,stepsize,tt,time]...
    = defineSimulationTimeFrames(printRequest);

%protons_per_ATPase = 3;     % three protons per ATP synthase (for output)
%mol2mmol = 1000;      

%%% Set default drug addition parameters (time and extent)
[rotenone, AA, oligo, CIV, FCCP, energy] = defineDefaultDrugCond(t_no_time);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Set what figures to plot
plot_ss = 0;    % Set to 1 to plot steady-state simulations
%%% Assign legends for output graphs 
[statevar_legend, output_legend] = getLegends; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate initial conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ctot0       = xpar(2);           % Cyt-c from Beard. Total IMS Cyt-c, Cred+Cox, molar
Qtot0       = xpar(3);           % total IMS Ubiquinol, Q+QH2, molar
ADTP_tot    = xpar(4);           % total Adenosine phosphates in cytosol, calculated from Evans_77

%%% Basal Respiratory State
% Set otherpar(4) (K_ADTP_dyn) = 0 to simulate no glycolytic ATP production (e.g. pyruvate buffer)
%otherpar(4) = 3; % Glycolytic capacity
state_fact  = 10/11; % 10/11 =) ss ATP:ADP = 10:1. 20/21 =) 20:1  3/4 =) 3:1 (State-3 mito)
ATP_e       = state_fact*ADTP_tot;  % ATP conc
ADP_e       = ADTP_tot-ATP_e;       % ADP conc
fprintf('Respiratory State (ATP:ADP): %0.1i\n', state_fact)
fprintf('Glycolytic Capacity (K_ADTP_dyn): %0.1i\n', otherpar(4))

xo_single_cell  = initial(ADP_e, ATP_e, Ctot0, Qtot0); % Initial2 contains updated ss values

%%% ODE options
options = odeset('RelTol',1e-5, 'AbsTol',1e-8, 'MaxStep',10e-1, ...
    'InitialStep',1e-1, 'MaxOrder',5, 'BDF','on');

% Define experiments to simulate (see defineExptsToSimulate.m)
num_expt = 1; % Initiaise value so baseline values only plotted once
expts = 5;    
% Set number of simulations (= number of cells in population)
numSims = 10;
printSim = 1;   % Output every simulation number to screen

for s = 1:expts
%if s == 1 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Steady-state Calculations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % oligo.check = 0;    % To ensure oligo doesn't get reset in ss calcs
  % Run sub-routine for set period to get steady-state 
  %[t0, y_ss]= ode15s(@sub_energetic,[t_prior t_start],xo_single_cell,options,xpar,otherpar,time,oligo,rotenone,AA,CIV,FCCP,energy);
  %fprintf('Steady-state calculated.\n')
  %fprintf('deltaPsi_m = %0.2i.\n',y_ss(end,19))
  %fprintf('ATP_c = %0.2i; ADP_c = %0.2i\n',y_ss(end,23),y_ss(end,24))
  %fprintf('ATP_x = %0.2i; ADP_x = %0.2i\n',y_ss(end,8),y_ss(end,9))
  %fprintf('NADH_x = %0.2i\n',y_ss(end,4))
%end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% ODE Calculations
  %%% Take last value of ss simulations as initial concentrations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Define experiment to simulate (based on input, s)
  printRequest = 1; %Set to 1 if you want to print info
  [rotenone,AA,oligo,CIV,FCCP,energy] = ...
      defineExptsToSimulate(s,t_no_time,printRequest);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set impairments to simulate
% Set parameters to 1 for 100% condition, and to <1 to
% simulate impairment
% Set v_otherpar4 = 0 to simulate pyruvate conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
v_xpar7 = 1;    % x_DH impairment: 0.7~99%, 0.55~97%, 0.4~90%, 0.32~80, 0.28~71%
v_xpar9 = 1;      % x_C1  
v_xpar10 = 1;     % x_C3: 5e-2
v_xpar13 = 1;     % x_C4 impairment: 5e-3~90%,3e-3~79%,2.1e-3~70%
v_xpar15 = 1;     % x_F1 impairment: 1e-5~70%; 2e-5~83%; 4.5e-5~92%
v_xpar21 = 1;     % x_Hle  
v_otherpar4 = 1;  % K_ADTP_dyn impairment: 0.9~93%, 0.7~79%, 0.6~70%
v_otherpar6 = 1;  %K_ADTP_cons

fprintf('\n********************')
if v_xpar7 < 1, fprintf('\nImpairing x_DH [xpar(7)=%0.2i]\n',v_xpar7), end
if v_xpar9 < 1, fprintf('\nImpairing x_C1 [xpar(9)=%0.2i]\n',v_xpar9), end
if v_xpar10 < 1, fprintf('\nImpairing x_C3 [xpar(10)=%0.2i]\n',v_xpar10), end
if v_xpar13 < 1, fprintf('\nImpairing x_C4 [xpar(13)=%0.2i]\n',v_xpar13), end
if v_xpar15 < 1, fprintf('\nImpairing x_F1 [xpar(15)=%0.2i]\n',v_xpar15), end
if v_xpar21 < 1, fprintf('\nImpairing x_Hle [xpar(21)=%0.2i]\n',v_xpar21), end
if v_otherpar4 < 1, fprintf('\nImpairing K_ADTP_dyn [otherpar(4)=%0.2i]\n',v_otherpar4), end
if v_otherpar6 < 1, fprintf('\nImpairing K_ADTP_cons [otherpar(6)=%0.2i]\n',v_otherpar6), end
fprintf('********************\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run multiple simulations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[stateVarAll,stateVarFCAll,stateVarAll1,...
    stateVarAll8,stateVarAll30,stateVarAll45,...
    stateVarAll70, OutputAll,OutputFCAll,OutputAll1,...
    OutputAll8,OutputAll30,OutputAll45,OutputAll70,...
    median_stateVarAll,median_stateVarFCAll,median_stateVarAll1,...
    median_stateVarAll8,median_stateVarAll30,median_stateVarAll45,...
    median_stateVarAll70,median_OutputAll,median_OutputFCAll,...
    median_OutputAll1,median_OutputAll8,median_OutputAll30,...
    median_OutputAll45,median_OutputAll70,...
    mean_stateVarAll,mean_stateVarFCAll,mean_stateVarAll1,...
    mean_stateVarAll8,mean_stateVarAll30,mean_stateVarAll45,...
    mean_stateVarAll70, mean_OutputAll,mean_OutputFCAll,mean_OutputAll1,...
    mean_OutputAll8,mean_OutputAll30,mean_OutputAll45,mean_OutputAll70] ...
    = populationRun(numSims, oligo,rotenone,AA,CIV,FCCP,energy,time,...
    v_xpar7,v_xpar9,v_xpar10,v_xpar13,v_xpar15,v_xpar21,...
    v_otherpar4,v_otherpar6,...
     t_prior,t_start,t_final,tt,stepsize,...
     options, xo_single_cell, plot_ss, statevar_legend, printSim);  
 
  % Save solutions to 3d matrices (so you can plot it similarly to
  % functions with multiple simulations
  median_stateVarAll_expts(:,:,s) = median_stateVarAll;
  median_stateVarFCAll_expts(:,:,s) = median_stateVarFCAll;
  median_OutputAll_expts(:,:,s) = median_OutputAll;
  median_OutputFCAll_expts(:,:,s) = median_OutputFCAll;
  stateVarAll1_expts(:,:,s) = stateVarAll1;
  OutputAll1_expts(:,:,s) = OutputAll1;
end   

%%%%%%%%%%%%%%%%%%%%
% Plot (median) baseline data vs literature values
if num_expt == 1
    % Only plot baseline values once
    num_expt = 2;
    type = 'population';
    plotBaselineModelVsLit(xpar,stateVarAll1(:,:,1),OutputAll1(:,:,1),type)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot time-course data vs Padova experiments
plotSingleRun_check_Timeseries(stepsize,tt,...
    median_OutputAll_expts, median_OutputFCAll_expts,...
    median_stateVarAll_expts, median_stateVarFCAll_expts)

title('**Median Values**')
end