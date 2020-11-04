function electrode_selection_threshold = compute_electrode_selection_threshold


clear;
clc;
username = char(java.lang.System.getProperty('user.name'));
toolpth_ = [ filesep 'Users' filesep username filesep 'Documents' filesep 'tools'];
pth_ = pwd;
assert(strcmp(pth_(end-3:end), 'code'), 'Please cd to folder');
data_pth_ = [pth_(1:end-4)];
addpath([toolpth_ filesep 'fieldtrip']);

ft_defaults;
sjs = [''; ''; '';''; ''; ''; ''; ''; '']; % codes removed for sharing



%% Loading RESULT
all_dat = {};
all_models = cell(size(sjs, 1),1);
all_mis = cell(size(sjs, 1),1);
for ss = 1 : size(sjs, 1)
    
    %subject code
    code_ = sjs(ss,:);
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    deriv_pth3 = [data_pth_ 'derivatives' filesep 'mi', filesep ...
        'sub-' code_ ];
    
    %load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max.mat']));
    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));
    try
        load(fullfile(deriv_pth3, ['MI_power_toHC.mat']));
    catch
        eeg_mi = [];
    end
    all_dat = cat(1,all_dat, model_eeg.effect);
    all_models{ss} = model_eeg;
    all_mis{ss} = eeg_mi;
end

%% compute the threshold: select electrodes based on GMM fitting

figure;

all_eff = cell2mat(all_dat);
GMModel = fitgmdist(all_eff, 2);

P = posterior(GMModel,all_eff);

[mu1, small_modelidx]  = min(GMModel.mu);
[mu2, big_modelidx]  = max(GMModel.mu);

selB = ((P(:, big_modelidx) ./ P(:,small_modelidx))>10)& (all_eff>0);

selA = ~selB;

y = pdf(GMModel,[-max(all_eff):0.00001:max(all_eff)]');


paramEsts= GMModel;
MU=[mu1;mu2];
SIGMA=cat(3,[paramEsts.Sigma(small_modelidx)],[paramEsts.Sigma(big_modelidx)]);
PPp=[paramEsts.PComponents(small_modelidx),paramEsts.PComponents(big_modelidx)];

objA = gmdistribution(MU(1),SIGMA(1),PPp(1));
xgridss=transpose(linspace(-max(all_eff),max(all_eff),1000));
tmp = pdf(objA,xgridss);

plot(xgridss,tmp./(sum(tmp)),'-', 'color',[30 83 145]./255, 'linewidth',2)
hold on;
objB = gmdistribution(MU(2),SIGMA(2),PPp(2));
xgridss=transpose(linspace(-max(all_eff),max(all_eff),1000));
tmp = pdf(objB,xgridss);
plot(xgridss,tmp./(sum(tmp)),'-','color',[198 0 27]./255, 'linewidth',2);

hold on;

plot([-max(all_eff):0.00001:max(all_eff)],y./sum(y), 'color',[39 40 56]./255,'linewidth',2);
xlim([-max(all_eff), max(all_eff)])
%ylim([0 max(y)+1000]);
line([min(all_eff(selB)),min(all_eff(selB))], [min(all_eff(selB)), 0.02], 'LineStyle', '--', 'LineWidth', 2, 'Color', [198 0 27]./255);


electrode_selection_threshold = min(all_eff(selB));

hold on
bp = boxplot(all_eff,'orientation', 'horizontal','positions', max(y./sum(y))+0.005,'colors', [39 40 56]./255, 'Widths',0.0025);
ylim([0 max(y./sum(y))+0.01]);
set(bp,'linewidth',2);
h=findobj(gca,'tag','Outliers');
set(h,'Visible', 'off')
hold on

plot(all_eff(all_eff<electrode_selection_threshold)...
    ,0.02*ones(size(...
    all_eff(all_eff<electrode_selection_threshold)))-0.0005.*rand(size( all_eff(all_eff<electrode_selection_threshold))),...
    'o','Color', [39 40 56]./255, 'MarkerSize',2, 'MarkerFaceColor', [39 40 56]./255);
hold on;
plot(all_eff(all_eff>=electrode_selection_threshold)...
    ,0.02*ones(size(...
    all_eff(all_eff>=electrode_selection_threshold)))-0.0005.*rand(size( all_eff(all_eff>=electrode_selection_threshold))),...
    'o','Color', [198 0 27]./255, 'MarkerSize',2, 'MarkerFaceColor', [198 0 27]./255);

set(gcf, 'color', 'w');
set(gca, 'ytick', []);
ylabel('density');
xlabel('value');




%%
end