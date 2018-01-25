function plotSingleRun_check_Timeseries(stepsize,tt,...
    OutputAll, OutputFCAll, stateVarAll, stateVarFCAll)


%%%%%%%%%%%%%%%%%%%%%
% Convert simulations to mol X/min/ug protein to align with exptal data
Output2Convert = OutputAll;
dim = 1; % dimension of Output2Convert variations
[OCR_converted,Output27Conv,Output5Conv,...
    J_ATP_total,prop_cyto,prop_mito,convert] = unitsConvert(Output2Convert,dim);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot time-course data vs Padova experiments
figure('Position',[100 100 1500, 800]);

%%%% OCR foldchange
subplot(2,4,1), hold on
%plot(tt,OutputFCAll(:,4,1)) % OCR during Rotenone treatment
%plot(tt,OutputFCAll(:,4,2)) % OCR during AA treatment
plot(tt(2:6/stepsize:end),OutputFCAll(2:6/stepsize:end,4,3)) % OCR during seahorse

% Obtain and plot experimenal data (Seahorse Glc FC)
[time, OCR_FC_S1, OCR_FC_S123, OCR_FC_S2, OCR_FC_S3,...
    OCR_norm_S1, OCR_norm_S123, OCR_norm_S2, OCR_norm_S3,...
    OCR_norm_IncEnDem] = exptData_OCR;

plot(time,OCR_FC_S1,'x')  
plot(time,OCR_FC_S123,'x') 
plot(time,OCR_FC_S2,'x')
plot(time,OCR_FC_S3,'x')
title('J-C4 (Padova vs sims)')
axis([0 inf -1 4])
ylabel('OCR (FC over baseline)')

% OCR raw (and converted) units
subplot(2,4,2), hold on
%plot(tt,OCR_converted(:,1)) % OCR during Rotenone treatment
%plot(tt,OutputFCAll(:,2)) % OCR during AA treatment
plot(tt(2:6/stepsize:end),OCR_converted(2:6/stepsize:end,3))
% Experimental data
plot(time,OCR_norm_S1,'x')
plot(time,OCR_norm_S123,'x')
plot(time,OCR_norm_S2,'x')
plot(time,OCR_norm_S3,'x')
ylabel('OCR (pmol O2/min/ug protein)')
title('OCR (- non-mito respiration)')
axis([0 inf -inf inf])

yyaxis right % activate right y-axis
plot(tt(2:6/stepsize:end),OutputAll(2:6/stepsize:end,4,3),'--') 
ylabel('J-C4 flux (mol O2/s/l mito)')
%legend('Sim. Glc OCR','OCR glcWT (mean) S1', 'S1+S2+S3', 'S2', 'S3', 'Sim. J-C4')

% OCR in increased energy demand conditions
subplot(2,4,3), hold on
plot(tt(2:6/stepsize:end),OCR_converted(2:6/stepsize:end,4))
% Rest color order of plot
ax = gca;
ax.ColorOrderIndex = 1;
plot(time,OCR_norm_IncEnDem,'x')
ylabel('OCR (pmol O2/min/ug protein)')
title('OCR IncEnDem (Picro - non-mito respiration)')
yyaxis right % activate right y-axis
plot(tt(2:6/stepsize:end),OutputAll(2:6/stepsize:end,4,4),'--',...
    'Color',[0.3010 0.7450 0.9330]) % Plot in light blue
ylabel('J-C4 flux (mol O2/s/l mito)')
%legend('Sim. Glc','OCR glc WT (mean)','Sim. J-C4 Glc')
% Set 2nd y-axis colour to match plotted line
ax = gca;
ax.YColor = [0.3010 0.7450 0.9330];

%%%%
%%%% TMRM, dPsi_m data
subplot(2,4,4),hold on
% Plot dPsi_m simulations (foldchange) for Rotenone and AntiA
plot(tt,stateVarFCAll(:,19,1))
plot(tt,stateVarFCAll(:,19,2))
%plot(tt,stateVarFCAll(:,19,3))

% Experimental TMRM data (Padova)
[time,ROF_Glc_WT_FC,AOF_Glc_WT_FC,ROF_Glc_WT_mV,AOF_Glc_WT_mV]...
    = exptData_TMRM;    % Obtain experimenal data
% Reset color order of plot
ax = gca;
ax.ColorOrderIndex = 1;
plot(time,ROF_Glc_WT_FC,'x') 
plot(time,AOF_Glc_WT_FC,'x')
ylabel('TMRM (\Delta\Psi_m) (FC over baseline)')

% activate right y-axis to plot mV values
yyaxis right 
% Simulations
plot(tt,stateVarAll(:,19,1),'--','Color',[0 0.447 0.741])
plot(tt,stateVarAll(:,19,2),'--','Color',[0.85 0.325 0.098])
% Experimental data
plot(time,ROF_Glc_WT_mV,'.','Color',[0 0.447 0.741]) 
plot(time,AOF_Glc_WT_mV,'.','Color',[0.85 0.325 0.098])
% Set 2nd y-axis colour to match plotted line
ax = gca;
ax.YColor = [0 0.447 0.741];
ylabel('\Delta\Psi_m (mV)')

%legend('Sim. ROF','Sim. AOF','ROFglc WT','AOFglc WT','raw Sim. ROF','raw Sim. AOF')
title('\Delta\Psi_m (Sims: FC- mV--, Expts: FCx mV_o)')

%%%%
%%%% NADH
subplot(2,4,5), hold on
plot(tt(1:50/stepsize),stateVarFCAll(1:50/stepsize,4,1))   % Rotenone
plot(tt(1:50/stepsize),stateVarFCAll(1:50/stepsize,4,3)) % Oligo+FCCP
%plot(tt(1:15/stepsize),stateVarFCAll(1:15/stepsize,4,4)) % IncEnDem
plot(tt(1:50/stepsize),stateVarFCAll(1:50/stepsize,4,5)) % FCCP, Rot

% Obtain and plot experimental data
[timeROF, timeOCR, timeIncEnDem, timeFR, NADH_Glc_ROF_FC,...
    NADH_Glc_OF_FC, NADH_Glc_IncEnDem_FC, NADH_Glc_FR_FC,...
    NADH_Glc_ROF_cal, NADH_Glc_OF_cal, NADH_Glc_IncEnDem_cal,...
    NADH_Glc_FR_cal] = exptData_NADH;

% Reset color order of plot
ax = gca;
ax.ColorOrderIndex = 1;
plot(timeROF,NADH_Glc_ROF_FC,'x')
plot(timeOCR,NADH_Glc_OF_FC,'x')
% Align expts with simulations (so only plot 3 mins before IncEnDem)
%plot(timeIncEnDem(10:end),NADH_Glc_IncEnDem_FC(10:end),'x')
plot(timeFR, NADH_Glc_FR_FC,'x')

% Plot calibrated foldchange values in same colours
ax = gca;
ax.ColorOrderIndex = 1;
plot(timeROF,NADH_Glc_ROF_cal,'.')
plot(timeOCR,NADH_Glc_OF_cal,'.')
%plot(timeIncEnDem(10:end),NADH_Glc_IncEnDem_cal(10:end),'.')
plot(timeFR, NADH_Glc_FR_cal,'.')
axis([0 50 0 2])
ylabel('NADH (FC over baseline)')

% Plot simulation raw values
yyaxis right % activate right y-axis
plot(tt(1:50/stepsize),stateVarAll(1:50/stepsize,4,1),'--',...
    'Color',[0 0.447 0.741])
plot(tt(1:50/stepsize),stateVarAll(1:50/stepsize,4,3),'--',...
    'Color',[0.85 0.325 0.098])
%plot(tt(1:15/stepsize),stateVarAll(1:15/stepsize,4,4),'--',...
%    'Color',[0.93 0.69 0.125])
plot(tt(1:50/stepsize),stateVarAll(1:50/stepsize,4,5),'--',...
    'Color',[0.93 0.69 0.125]);  %[0.49 0.18 0.56])
% Set 2nd y-axis colour to match plotted line
ax = gca;
ax.YColor = [0 0.447 0.741];
ylabel('[NADH] (M)')

title('NADH (ROF,OF,FR)')
%legend('Sim. NADH Glc ROF','NADH glcWT-ROF','Raw Sim. ROF',...
%    'Location','southwest')

%%%%% J_ATP_cyto, J_ATP_mito, J_ATP_total 
subplot(2,4,7), hold on
plot(tt(2:4/stepsize:40/stepsize),Output27Conv(2:4/stepsize:40/stepsize,3),'b') % J_ATP_cyto
plot(tt(2:4/stepsize:40/stepsize),Output5Conv(2:4/stepsize:40/stepsize,3),'r')  % J_ATP_mito
plot(tt(2:4/stepsize:40/stepsize),J_ATP_total(2:4/stepsize:40/stepsize,3),'g') % J_ATP_total
% Obtain exptal data
[time, J_ATP_cyto, J_ATP_mito, J_ATP_tot,...
    J_ATP_cyto_prop, J_ATP_mito_prop] = exptData_ATP; 
plot(time,J_ATP_cyto,'bx')
plot(time,J_ATP_mito,'rx')
plot(time,J_ATP_tot,'gx')
ylabel('J-ATP-prod (pmol ATP/min/ug protein)')
yyaxis right % activate right y-axis
plot(tt(2:4/stepsize:40/stepsize),OutputAll(2:4/stepsize:40/stepsize,27,3),'b--')
plot(tt(2:4/stepsize:40/stepsize),OutputAll(2:4/stepsize:40/stepsize,5,3),'r--')
ATP_total = OutputAll(:,27,3) + OutputAll(:,5,3);
plot(tt(2:4/stepsize:40/stepsize),ATP_total(2:4/stepsize:40/stepsize),'g--')
ylabel('Sim J-ATP-prod (mol ATP/s/l mito')
%legend('Sim cyto','mito','total','Padova cyto','mito','total','Raw sim cyto','mito','total')

%%% Plot proportions of cyto and mito ATP production
subplot(2,4,8), hold on
% Simulations
plot(tt(2:4/stepsize:40/stepsize),prop_cyto(2:4/stepsize:40/stepsize,3))
plot(tt(2:4/stepsize:40/stepsize),prop_mito(2:4/stepsize:40/stepsize,3))
% Experimental data
plot(time,J_ATP_cyto_prop,'bx')
plot(time,J_ATP_mito_prop,'rx')
ylabel('% contribution to total ATP prod')
legend('Sim. % cyto','% mito','Expt % cyto','% mito','Location','east')
title('% contrib. to total ATP (Oligo)')

end