function Beard_NC_singleRun_check
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Code adapted from crosst.m of Huber 2011 (MSB) and 
%%% regulation.m of Huber 2012 (MBS). 
%%% This function runs one simulation of basic functions
%%% Plots important steady-state parameter values (and literature values)
%%% and timecourses
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
numSims = 5;    

for s = 1:numSims
if s == 1 
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Steady-state Calculations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  oligo.check = 0;    % To ensure oligo doesn't get reset in ss calcs
  % Run sub-routine for set period to get steady-state 
  [t0, y_ss]= ode15s(@sub_energetic,[t_prior t_start],xo_single_cell,options,xpar,otherpar,time,oligo,rotenone,AA,CIV,FCCP,energy);
  fprintf('Steady-state calculated.\n')
  fprintf('deltaPsi_m = %0.2i.\n',y_ss(end,19))
  fprintf('ATP_c = %0.2i; ADP_c = %0.2i\n',y_ss(end,23),y_ss(end,24))
  fprintf('ATP_x = %0.2i; ADP_x = %0.2i\n',y_ss(end,8),y_ss(end,9))
  fprintf('NADH_x = %0.2i\n',y_ss(end,4))
end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% ODE Calculations
  %%% Take last value of ss simulations as initial concentrations
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Define experiment to simulate (based on input, s)
  printRequest = 1; %Set to 1 if you want to print info
  [rotenone,AA,oligo,CIV,FCCP,energy] = ...
      defineExptsToSimulate(s,t_no_time,printRequest);

  oligo.check     = 1;        % Oligo trigger (stores J_F1 value prior to Oligo)
  [t, y] = ode15s(@sub_energetic,[t_start t_final],y_ss(end,:),options,xpar,otherpar,time,oligo,rotenone,AA,CIV,FCCP,energy);

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Calculate outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Output = zeros(length(t),size(C,2));
  for i = 1:length(t)
    t_out = t(i); 
    y_out = y(i,:)'; 
    [f] = sub_energetic(t_out,y_out,xpar,otherpar,time,oligo,rotenone,AA,CIV,FCCP,energy);
    Output(i,:) = C;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Calculate fold changes of state variables and outputs
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Convert all solutions to set number of datapoints
  % Careful here because this samples the solution - may miss very rapid
  % changes
  for i = 1:size(y,2)    
    stateVarSpline(:,i) = spline(t, y(:,i), tt);    
  end
  for i = 1:size((Output),2)    
    OutputSpline(:,i) = spline(t, Output(:,i), tt);    
  end

  % Calculate fold change
  for i = 1:size(y,2)
    stateVarFC(:,i) = (stateVarSpline(:,i)./stateVarSpline(2,i)); 
    % Baseline is second value of solution (to ignore initial fluctuations)
  end
  for i = 1:size(Output,2)
    OutputFC(:,i) = (OutputSpline(:,i)./OutputSpline(2,i)); 
  end
  
  % Save solutions into 3d matrices (so you can plot it similarly to
  % functions with multiple simulations)
  stateVarAll(:,:,s) = stateVarSpline;
  stateVarFCAll(:,:,s) = stateVarFC;
  OutputAll(:,:,s) = OutputSpline;
  OutputFCAll(:,:,s) = OutputFC;
    
end   

%%%%%%%%%%%%%%%%%%%%
% Plot baseline data vs literature
type = 'singleSim';
plotBaselineModelVsLit(xpar,stateVarAll,OutputAll,type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot time-course data vs Padova experiments
plotSingleRun_check_Timeseries(stepsize,tt,...
    OutputAll, OutputFCAll, stateVarAll, stateVarFCAll)

end