function plotBaselineModelVsLit(xpar,stateVarAll,OutputAll,type)
% Plot baseline values next to literature values
% Literature values listed in "ICs and Ks from lit.xls"

figure('Position',[100 100 1500, 800],'Name', 'Baseline Values');

%%%
% dPsim
subplot(3,5,1)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,19,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),stateVarAll(:,19,1),'filled','k')
    hold on
    y = nanmedian(stateVarAll(:,19,1));
    z = nanmean(stateVarAll(:,19,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on, ylabel('Mito. memb pot \Delta\Psi_m (mV)')
% Gerencser 2012
scatter(1.2,139,'filled','b')
plot([1.15 1.25],[134 134],'b'), plot([1.15 1.25],[144 144],'b') 
% Ward 2000, Nicholls & Budd 2000, Nicholls 1980
scatter(1.3,150,'filled','b')
% Surin 2013
plot([1.35 1.45],[150 150],'b'), plot([1.35 1.45],[180 180],'b')
axis([0.9 inf -inf inf])

%%%
% Cytosolic [ATP]
subplot(3,5,2)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,23,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),stateVarAll(:,23,1),'filled','k')
    hold on
    y = nanmedian(stateVarAll(:,23,1));
    z = nanmean(stateVarAll(:,23,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
ylabel('Cyto [ATP] (M)'), hold on
scatter(1.2,2.68e-3,'filled','b') % Weisova
scatter(1.3,2.59e-3,'filled','b') % Lehninger
scatter(1.4,2.2e-3,'filled','m') % Aubert model
scatter(1.45,2.2e-3,'filled','m') % DiNuzzo model
scatter(1.5,2.25e-3,'filled','m') % Cloutier model
%plot([0.9 1.1], [2.2e-3 2.2e-3], 'b')
%plot([0.9 1.1], [2.7e-3 2.7e-3], 'b')
axis([0.9 inf -inf inf])

%%%
% Cytosolic [ADP]
subplot(3,5,3)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,24,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),stateVarAll(:,24,1),'filled','k')
    hold on
    y = nanmedian(stateVarAll(:,24,1));
    z = nanmean(stateVarAll(:,24,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on
scatter(1.3,730e-6,'filled','b') % Lehninger
ylabel('Cyto [ADP] (M) [ATP:ADP 10:1 Hardie]'), hold on
axis([0.9 inf -inf inf])

%%%
% Cytosolic ATP:ADP
subplot(3,5,4)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,23,1)/stateVarAll(2,24,1),'filled','k')
    case 'population'
    ratio = stateVarAll(:,23,1)./stateVarAll(:,24,1);
    scatter(ones(size(stateVarAll,1),1),ratio,'filled','k')
    hold on
    y = nanmedian(ratio);
    z = nanmean(ratio);
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on, scatter(1.2,10,'b','filled')    % Hardie, Nicholls 2000
ylabel('Cyto ATP/ADP')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial ATP
subplot(3,5,5)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,8,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),stateVarAll(:,8,1),'filled','k')
    hold on
    y = nanmedian(stateVarAll(:,8,1));
    z = nanmean(stateVarAll(:,8,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on
scatter(1.2,2.5e-3,'r','filled') % Yoshida PC12 (they measured ATP_c as 3.9 mM)
ylabel('Mito [ATP] (M) (*PC12))')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial [ADP]
subplot(3,5,6)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,9,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),stateVarAll(:,9,1),'filled','k')
    hold on
    y = nanmedian(stateVarAll(:,9,1));
    z = nanmean(stateVarAll(:,9,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on, ylabel('Mito [ADP] (M)')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial ATP:ADP
subplot(3,5,7)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,8,1)/stateVarAll(2,9,1),'filled','k'), hold on
    case 'population'
    ratio = stateVarAll(:,8,1)./stateVarAll(:,9,1);
    scatter(ones(size(stateVarAll,1),1),ratio,'filled','k'), hold on
    y = nanmedian(ratio);
    z = nanmean(ratio);
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
ylabel('Mito ATP/ADP')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial [NADH]
subplot(3,5,8)
switch type
    case 'singleSim'
        scatter(1,stateVarAll(2,4,1),'filled','k')
        case 'population'
        scatter(ones(size(stateVarAll,1),1),stateVarAll(:,4,1),'filled','k')
        hold on
        y = nanmedian(stateVarAll(:,4,1));
        z = nanmean(stateVarAll(:,4,1));
        plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on
plot([1.15 1.25], [25e-6 25e-6], 'm'), plot([1.15 1.25], [400e-6 400e-6], 'r')
ylabel('Mito [NADH] (M)')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial [NAD]
subplot(3,5,9)
switch type
    case 'singleSim'
        scatter(1,xpar(1)-stateVarAll(2,4,1),'filled','k')
        case 'population'
        diff = xpar(1)-stateVarAll(:,4,1);
        scatter(ones(size(stateVarAll,1),1),diff,'filled','k')
        hold on
        y = nanmedian(diff);
        z = nanmean(diff);
        plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on
plot([1.15 1.25],[250e-6 250e-6],'r'), plot([1.2 1.2],[250e-6 300e-6],'r')
ylabel('Mito [NAD] (M)')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial NAD:NADH
subplot(3,5,10)
switch type
    case 'singleSim'
    scatter(1,(xpar(1)-stateVarAll(2,4,1))/stateVarAll(2,4,1),'filled','k')
    case 'population'
    ratio = (xpar(1)-stateVarAll(:,4,1))./stateVarAll(:,4,1);
    scatter(ones(size(stateVarAll,1),1),ratio,'filled','k')
    hold on
    y = nanmedian(ratio);
    z = nanmean(ratio);
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
ylabel('Mito NAD:NADH')
hold on, scatter(1.2,9,'b','filled') % Nicholls & Budd 2000
plot([1.25 1.35],[7 7],'r'), plot([1.25 1.35],[8 8],'r') % Lots
plot([1.38 1.45],[1 1],'r'), plot([1.38 1.45],[8 8],'r') % Bilan 2014
plot([1.45 1.55],[4 4],'r'), plot([1.45 1.55],[10 10],'r') % link in .xls
plot([1.55 1.65],[2 2],'r'), plot([1.55 1.65],[16 16],'r') % Zhao 2014
axis([0.9 inf -inf inf])

%%%
% Mitochondrial pH [H+]
subplot(3,5,11)
switch type
    case 'singleSim'
    scatter(1,OutputAll(2,26,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),OutputAll(:,26,1),'filled','k')
    hold on
    y = nanmedian(OutputAll(:,26,1));
    z = nanmean(OutputAll(:,26,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on, scatter(1.2,8,'b','filled') % Bolshakov 2007
plot([1.25 1.35],[7.5 7.5],'r'),plot([1.25 1.35],[8 8],'r') % Poburko 2011
ylabel('Mito pH')
axis([0.9 inf -inf inf])

%%%
% Mitochondrial O2
subplot(3,5,12)
switch type
    case 'singleSim'
    scatter(1,stateVarAll(2,7,1),'filled','k')
    case 'population'
    scatter(ones(size(stateVarAll,1),1),stateVarAll(:,7,1),'filled','k')
    hold on
    y = nanmedian(stateVarAll(:,7,1));
    z = nanmean(stateVarAll(:,7,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
hold on, ylabel('Mitochondrial O_2 (M)')
scatter(1.2,26e-6,'m','filled') % Aubert model
scatter(1.3,175e-6,'b','filled') % Dmitriev 2015
plot([1.35 1.45],[3e-6 3e-6],'r')
plot([1.35 1.45],[30e-6 30e-6],'r') % Murphy 09
axis([0.9 inf -inf inf])

%%%%%%%%%%%
% Convert simulations to mol X/min/ug protein to align with exptal data
Output2Convert = OutputAll;
dim = 1; % dimension of Output2Convert variations
[OCR_converted,Output27Conv,Output5Conv,...
    J_ATP_total,prop_cyto,prop_mito] = unitsConvert(Output2Convert,dim);

% Converted basal OCR
subplot(3,5,13)
switch type
    case 'singleSim'
    scatter(1, OCR_converted(2,1,1), 'filled','k'), hold on
    case 'population'
    scatter(ones(size(OCR_converted,1),1), OCR_converted(:,1,1), 'filled','k')
    hold on
    y = nanmedian(OCR_converted(:,1,1));
    z = nanmean(OCR_converted(:,1,1));
    plot([0.95 1.05],[y y],'k'), plot([0.95 1.05],[z z],'k--')
end
% Different experimental series (Pierre)
scatter(1.2,5.1,'b','filled'), scatter(1.2,7.5,'b','filled')
scatter(1.2,9.5,'b','filled')
ylabel('OCR basal (mol O2/min/ug protein)')
axis([0.9 inf -inf inf])

%%%
% Converted cyto, mito and total ATP production
subplot(3,5,14), hold on
switch type
    case 'singleSim'
    scatter(1, Output27Conv(2,1,1), 'b', 'filled')
    scatter(1, Output5Conv(2,1,1), 'r', 'filled')
    scatter(1, J_ATP_total(2,1,1), 'y', 'filled')
    case 'population'
    scatter(ones(size(Output27Conv,1),1), Output27Conv(:,1,1), 'b', 'filled')
    scatter(ones(size(Output5Conv,1),1), Output5Conv(:,1,1), 'r', 'filled')
    scatter(ones(size(J_ATP_total,1),1), J_ATP_total(:,1,1), 'y', 'filled')
end
plot([0.9 1.1], [24 24], 'r'), plot([0.9 1.1], [15.5 15.5], 'b')
plot([0.9 1.1], [43.1 43.1], 'y')
ylabel('J-ATP (mol ATP/min/ug protein)')
legend('cyto','mito','total')
axis([0.9 inf -inf inf])

%%%
% Contribution of cyto and mito to total ATP production 
subplot(3,5,15), hold on
switch type
    case 'singleSim'
    scatter(1, prop_cyto(2,1,1), 'filled')
    scatter(1, prop_mito(2,1,1), 'filled')
    case 'population'
    scatter(ones(size(prop_cyto,1),1), prop_cyto(:,1,1), 'filled', 'b')
    scatter(ones(size(prop_cyto,1),1), prop_mito(:,1,1), 'filled', 'r')
    yCyto = nanmedian(prop_cyto(:,1,1));
    zCyto = nanmean(prop_cyto(:,1,1));
    yMito = nanmedian(prop_mito(:,1,1));
    zMito = nanmean(prop_mito(:,1,1));
    plot([0.95 1.05],[yCyto yCyto],'b'), plot([0.95 1.05],[zCyto zCyto],'b--')
    plot([0.95 1.05],[yMito yMito],'b'), plot([0.95 1.05],[zMito zMito],'b--')
end
% Exptal measurements
plot([0.9 1.1], [39 39], 'b'), plot([0.9 1.1], [61 61], 'r')
ylabel('Contribution to total ATP prod.')
legend('cyto','mito')
axis([0.9 inf -inf inf])
end