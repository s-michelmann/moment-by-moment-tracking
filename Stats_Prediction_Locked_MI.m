%%
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

% electrode_selection_threshold
electrode_selection_threshold = compute_electrode_selection_threshold;

%% load the data

epos_all = [];
data_all1 = [];
data_all2 = [];

data_all1cross = [];
data_all2cross = [];

drop_delta_all1 = [];
drop_delta_all2 = [];
drop_theta_all1 = [];
drop_theta_all2 = [];
drop_beta_all1 = [];
drop_beta_all2 = [];
drop_alpha_all1 = [];
drop_alpha_all2 = [];
drop_gamma_low_all1 = [];
drop_gamma_low_all2 = [];
drop_gamma_high_all1 = [];
drop_gamma_high_all2 = [];
drop_raw_all1 = [];
drop_raw_all2 = [];

single_delta_all1 = [];
single_delta_all2 = [];
single_theta_all1 = [];
single_theta_all2 = [];
single_beta_all1 = [];
single_beta_all2 = [];
single_alpha_all1 = [];
single_alpha_all2 = [];
single_gamma_low_all1 = [];
single_gamma_low_all2 = [];
single_gamma_high_all1 = [];
single_gamma_high_all2 = [];
single_raw_all1 = [];
single_raw_all2 = [];

label_all = [];
hc_idices = [];
eff_indices = [];
hc_labels = {};
n_hc_elecs = [];
codes_all = {};
for ss = 1 : size(sjs, 1)
    % get sj code
    code_ = sjs(ss,:);

    % find path names
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    deriv_pth3 = [data_pth_ 'derivatives' filesep 'mi', filesep ...
        'sub-' code_ ];
    % load the EEG
    try   
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksFinal1s_slidingWindow.mat']);       
    catch
         disp(['no data for sj ' code_ ])
        continue;
    end

    
    % load the GC results
    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));
    fclose all;
    fid = fopen([deriv_pth filesep 'HC_channels.tsv'], 'r');
    txt = textscan(fid, '%s', 'delimiter', '\t', 'Headerlines', 0);
    fclose all   ;
    hc_tmp = zeros(size(model_eeg.label));
    tmp_n = 0;
    % if there are no channels continue with the next subject
    if ~isempty(txt{1})
        hc_labels = cat(1, hc_labels, txt{:});
        hc_idx = find( cell2mat(cellfun(@(x) ismember(x, txt{:}), model_eeg.label, 'Un', 0)));
        hc_tmp(hc_idx) = 1;
        tmp_n = numel(hc_idx);
    end
    hc_idices = cat(1, hc_idices, hc_tmp);
    
    
    % keep the indices of the electrodess that show an effect above thresh
    eff_idx = (cell2mat(model_eeg.effect) >= electrode_selection_threshold);
    if ~any(eff_idx) 
        disp(['no effect for sj ' code_ ])
        continue;
        
    end    
    eff_indices = cat(1, eff_indices, eff_idx);
    
    
    % keep the indices of the channels
    n_hc_elecs = [n_hc_elecs; tmp_n];
    label_all = [label_all; model_eeg.label];
    load(fullfile(deriv_pth, 'elec_orig.mat'));
    epos = getEposFromElec(elec_orig, model_eeg.label);

    epos_all = [epos_all;epos];

    data_all1 = [data_all1; eeg_mi.trial{1}];
    data_all2 = [data_all2; eeg_mi.trial{2}];
    codes_all = [codes_all; code_]
    
    try
        load([deriv_pth3 filesep 'CrossLaggedMIaround_predictionPeaks1s_slidingWindow.mat']);
        data_all1cross = [data_all1cross; eeg_mi.trial{1}];
        data_all2cross = [data_all2cross; eeg_mi.trial{2}];
    catch
        disp(ss);
        continue;
    end
    
    try
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDropdelta.mat']);

        drop_delta_all1 = [drop_delta_all1 ; eeg_mi.trial{1}];
        drop_delta_all2 = [drop_delta_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDroptheta.mat']);
        drop_theta_all1 = [drop_theta_all1 ; eeg_mi.trial{1}];
        drop_theta_all2 = [drop_theta_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDropbeta.mat']);
        drop_beta_all1 = [drop_beta_all1 ; eeg_mi.trial{1}];
        drop_beta_all2 = [drop_beta_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDropalpha.mat']);
        drop_alpha_all1 = [drop_alpha_all1 ; eeg_mi.trial{1}];
        drop_alpha_all2 = [drop_alpha_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDroplow_gamma.mat']);
        drop_gamma_low_all1 = [drop_gamma_low_all1 ; eeg_mi.trial{1}];
        drop_gamma_low_all2 = [drop_gamma_low_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDrophigh_gamma.mat']);
        drop_gamma_high_all1 = [drop_gamma_high_all1 ; eeg_mi.trial{1}];
        drop_gamma_high_all2 = [drop_gamma_high_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksDropraw.mat']);
        drop_raw_all1 = [drop_raw_all1 ; eeg_mi.trial{1}];
        drop_raw_all2 = [drop_raw_all2 ; eeg_mi.trial{2}];
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSingledelta.mat']);
        single_delta_all1 = [single_delta_all1 ; eeg_mi.trial{1}];
        single_delta_all2 = [single_delta_all2 ; eeg_mi.trial{2}];
        
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSingletheta.mat']);
        single_theta_all1 = [single_theta_all1 ; eeg_mi.trial{1}];
        single_theta_all2 = [single_theta_all2 ; eeg_mi.trial{2}];
        
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSinglebeta.mat']);
        single_beta_all1 = [single_beta_all1 ; eeg_mi.trial{1}];
        single_beta_all2 = [single_beta_all2 ; eeg_mi.trial{2}];
        
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSinglealpha.mat']);
        single_alpha_all1 = [single_alpha_all1 ; eeg_mi.trial{1}];
        single_alpha_all2 = [single_alpha_all2 ; eeg_mi.trial{2}];
        
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSinglelow_gamma.mat']);
        single_gamma_low_all1 = [single_gamma_low_all1 ; eeg_mi.trial{1}];
        single_gamma_low_all2 = [single_gamma_low_all2 ; eeg_mi.trial{2}];
        
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSinglehigh_gamma.mat']);
        single_gamma_high_all1 = [single_gamma_high_all1 ; eeg_mi.trial{1}];
        single_gamma_high_all2 = [single_gamma_high_all2 ; eeg_mi.trial{2}];
        
        
        load([deriv_pth3 filesep 'LaggedMIaround_predictionPeaksSingleraw.mat']);
        single_raw_all1 = [single_raw_all1 ; eeg_mi.trial{1}];
        single_raw_all2 = [single_raw_all2 ; eeg_mi.trial{2}];

 
    catch
        disp(ss);
        continue;
    end

end
%% plot network in color

yeo_colors = [0	 0	0 ;...
    120	18	134; ...
    255	0	0	;...
    70	130	180;...
    42	204	164;...
    74	155	60;...
    0	118	14;...
    196	58	250;...
    255	152	213;...
    200	248	164;...
    122	135	50;...
    119	140	176;...
    230	148	34;...
    135	50	74;...
    12	48	255;...
    0	0	130;...
    255	255	0;...
    205	62	78]./255;

[label_nr] = get_YEO_labels(epos_all);

label_nr(find(eff_indices)) = -1;
label_nr(find(hc_idices)) = 18;
% figure;
plot_ecog(ones(size(find(label_nr ==1))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos_all(find(label_nr ==1),:),...
    [], 0.2, [-90 0], 22, 1, 2, repmat(yeo_colors(1+1,:), [256 1]));

tst = find(label_nr ==1);
view([-90 0]);
hold on;
for ll = 2 : max(label_nr)-1
    tst = [tst; find(label_nr ==ll)];
    plot_ecog(ones(size(find(label_nr ==ll))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
        ,epos_all(find(label_nr ==ll),:),...
        [], 0.2, [-90 0], 22, 0, 2, repmat(yeo_colors(ll+1,:), [256 1]));   
end
set(gcf, 'color', 'white');
% 
% figure;
% plot_ecog(ones(size(find(label_nr))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
%         ,epos_all(find(label_nr),:),...
%         [], 0.2, [-90 0], 22, 0, 2, repmat([1 0 0 ], [256 1]));   

%% plot the hippocampus
n_ = numel(find(hc_idices));
yl_= 0.045;
figure;
subplot(141);
plot_ecog(ones(size(find(hc_idices))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos_all(find(hc_idices),:),...
    [], 0.2, [0 90], 22, 1, 2, repmat([1 0 0], [256 1]));

subplot(142);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all1(find(hc_idices),:)) , ...
    nanstd(data_all1(find(hc_idices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([0 yl_]);
title('run1')

subplot(143);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(hc_idices),:)) , ...
    nanstd(data_all2(find(hc_idices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([0 yl_]);
title('run2');

[H, p , CI, stats] = ttest(data_all2(find(hc_idices),:), ...
    data_all1(find(hc_idices),:));

subplot(144);
plot(eeg_mi.time{1}.*10, stats.tstat); hold on
plot(eeg_mi.time{1}.*10, ones(size(stats.tstat)).*tinv(0.95, n_-1), 'r');
plot(eeg_mi.time{1}.*10, ones(size(stats.tstat)).*tinv(0.05, n_-1), 'r')
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([0 0.045]);
title('run2');
xlim([-1500 1500]);
ylim([-8 8]);
set(gcf, 'Position',  [600, 600, 1400, 400]);
p(301) = 1;
pvec = p(150:451);
p_thrsh = fdr(pvec, 0.05)
%% pretty plot plus stats on the hippocampus


n_ = numel(find(hc_idices));
yl_= 0.045;
[H, p , CI, stats] = ttest(data_all2(find(hc_idices),:), ...
    data_all1(find(hc_idices),:));
plim_ = fdr(p, 0.05);

figure;
subplot(131);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all1(find(hc_idices),:)) , ...
    nanstd(data_all1(find(hc_idices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([0 yl_]);

title('run1')
subplot(132);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(hc_idices),:)) , ...
    nanstd(data_all2(find(hc_idices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([0 yl_]);
title('run2');

subplot(133);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(hc_idices),:)- data_all1(find(hc_idices),:)) , ...
    nanstd(data_all2(find(hc_idices),:)-data_all1(find(hc_idices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Difference in Mutual Information');
xlim([-1500 1500]);
% ylim([0 0.045]);
title('run2 - run1');

set(gcf, 'Position',  [600, 600, 1400, 400]);

hold on;
tmp = nan(size(eeg_mi.time{1}));
tmp(p<plim_) = 0 ;
plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',8);

%% plot other


n_ = numel(find(~hc_idices&~eff_indices));
yl_= 0.045;
[H, p , CI, stats] = ttest(data_all2(find(~hc_idices&~eff_indices),:), ...
    data_all1(find(~hc_idices&~eff_indices),:));
plim_ = fdr(p, 0.001);

figure;
subplot(131);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all1(find(~hc_idices&~eff_indices),:)) , ...
    nanstd(data_all1(find(~hc_idices&~eff_indices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-3000 3000]);
ylim([0 yl_]);

title('run1')
subplot(132);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(~hc_idices&~eff_indices),:)) , ...
    nanstd(data_all2(find(~hc_idices&~eff_indices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-3000 3000]);
ylim([0 yl_]);
title('run2');

subplot(133);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(~hc_idices&~eff_indices),:)- data_all1(find(~hc_idices&~eff_indices),:)) , ...
    nanstd(data_all2(find(~hc_idices&~eff_indices),:)-data_all1(find(~hc_idices&~eff_indices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Difference in Mutual Information');
xlim([-3000 3000]);
% ylim([0 0.045]);
title('run2 - run1');

set(gcf, 'Position',  [600, 600, 1400, 400]);

hold on;
tmp = nan(size(eeg_mi.time{1}));
tmp(p<plim_) = 0 ;
plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',8);

%% plot ROIS

% exclude hc_indices
other_labels = label_nr(find(~hc_idices&~eff_indices));
other_epos = epos_all(find(~hc_idices&~eff_indices),:);
hc_epos = epos_all(find(hc_idices&~eff_indices),:);
data_other1 = data_all1(find(~hc_idices&~eff_indices),:);
data_other2 = data_all2(find(~hc_idices&~eff_indices),:);


for ll = 1 : 17
    n_ = numel(find(other_labels ==ll));
    dt1 = nanmean(data_other1(find(other_labels ==ll),:));
    dt2 = nanmean(data_other2(find(other_labels ==ll),:));
    yl_= max([dt1, dt2]);
    figure;
    subplot(141);
    plot_ecog(ones(size(find(other_labels ==ll))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
        ,other_epos(find(other_labels ==ll),:),...
        [], 0.2, [0 90], 22, 1, 2, repmat(yeo_colors(ll+1,:), [256 1]));
   
    subplot(142);
    shadedErrorBar(eeg_mi.time{1}.*10, ...
        nanmean(data_other1(find(other_labels ==ll),:)) , ...
        nanstd(data_other1(find(other_labels ==ll),:))./sqrt(n_));
    xlabel('HC shift relative to peak (ms)');
    ylabel('Mutual Information');
    xlim([-3000 3000]);
    ylim([0 yl_]);
    title('run1')
    
    subplot(143);
       shadedErrorBar(eeg_mi.time{1}.*10, ...
        nanmean(data_other2(find(other_labels ==ll),:)) , ...
        nanstd(data_other2(find(other_labels ==ll),:))./sqrt(n_));
    xlabel('HC shift relative to peak (ms)');
    ylabel('Mutual Information');
    xlim([-3000 3000]);
    ylim([0 yl_]);
    title('run2');
    
    [H, p , CI, stats] = ttest(data_other2(find(other_labels ==ll),:), ...
        data_other1(find(other_labels ==ll),:));
    
    subplot(144);
    plot(eeg_mi.time{1}.*10, stats.tstat); hold on
    plot(eeg_mi.time{1}.*10, ones(size(stats.tstat)).*tinv(0.95, n_-1), 'r');
    plot(eeg_mi.time{1}.*10, ones(size(stats.tstat)).*tinv(0.05, n_-1), 'r')
    xlabel('HC shift relative to peak (ms)');
    ylabel('Mutual Information');
    xlim([-3000 3000]);
    ylim([0 0.045]);
    title('run2');
    xlim([-3000 3000]);
    ylim([-8 8]);
    set(gcf, 'Position',  [600, 600, 1400, 400])
end

%%
n_ = size(data_all1(find(~hc_idices&~eff_indices)),1);
figure;
subplot(121);
title('run1')
shadedErrorBar(eeg_mi.time{1}.*10, nanmean(data_all1(find(~hc_idices&~eff_indices),:)), ...
    nanstd(data_all1(find(~hc_idices&~eff_indices),:))./sqrt(n_));
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-3000 3000]);
ylim([0 0.045]);
subplot(122);

title('run2')
shadedErrorBar(eeg_mi.time{2}.*10, nanmean(data_all2(find(~hc_idices&~eff_indices),:)), ...
    nanstd(data_all2(find(~hc_idices&~eff_indices),:))./sqrt(n_));
ylim([0 0.045]);
xlim([-3000 3000]);
xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
%% T test the differece vs. other

[H, p , CI, stats] = ttest2(data_all2(find(hc_idices),:) - ...
    data_all1(find(hc_idices),:), data_all2(find(~hc_idices&~eff_indices),:) - ...
    data_all1(find(~hc_idices&~eff_indices),:));


plot(eeg_mi.time{1}.*10, stats.tstat); hold on
plot(eeg_mi.time{1}.*10, ones(size(stats.tstat)).*tinv(0.95, n_-1), 'r');
plot(eeg_mi.time{1}.*10, ones(size(stats.tstat)).*tinv(0.05, n_-1), 'r')
xlabel('HC shift relative to peak (ms)');
ylabel('t-value');


plim_ = fdr(p, 0.05);
tmp = nan(size(eeg_mi.time{1}));
tmp(p<plim_) = 0 ;
plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',8);

title('difference in HC vs. in Other (independent sample t-test')
% %%
%
%
% %%
%
% title('difference')
% shadedErrorBar(eeg_mi.time{2}.*10, nanmean(data_all2-data_all1), ...
%     nanstd(data_all2-data_all1)./sqrt(30));
% hold on;
% line([-3000 3000],[0 0], 'color', [0.7 0.5 0.5],'LineStyle','--', 'LineWidth',2);
% xlim([-3000 3000]);
% % xlim([-1000 1000]);
% xlabel('HC shift relative to peak (ms)');
% ylabel('Mutual Information');
%
% [H, p , CI, stats] = ttest(data_all2, data_all1);
%% now load in the peak-shuffled data
yeo_colors(19,:) = [1 0 0];
for roi_ = 18 %:18; %18 is HC
    
    data_all1rand = nan(numel(find((label_nr == roi_))), 601, 1000);
    data_all2rand = nan(numel(find((label_nr == roi_))), 601, 1000);
    
    % Hippocampus first!
    
    for rr = 1 : 1000
        rr
        data_all1randtmp = [];
        data_all2randtmp = [];
        for ss = 1 : size(sjs, 1)
            % get sj code
            code_ = sjs(ss,:);
            % find path names
            deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
                'sub-' code_ ];
            deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
                'sub-' code_ ];
            deriv_pth3 = [data_pth_ 'derivatives' filesep 'mi', filesep ...
                'sub-' code_ ];
            deriv_pth4 = [data_pth_ 'derivatives' filesep 'mi_shuffle', filesep ...
                'sub-' code_ ];
            if ~exist([deriv_pth4 filesep 'RandomLaggedMIaround_predictionPeaks1.mat']); continue; end
            
       
            % load the EEG
            try
                
                load([deriv_pth4 filesep 'RandomLaggedMIaround_predictionPeaks'  num2str(rr) '.mat']);
            catch
                disp(ss);
                continue;
            end
            
            %load([deriv_pth filesep 'eeg_manualica_notch.mat']);
            
            data_all1randtmp = [data_all1randtmp; eeg_mi.trial{1}];
            data_all2randtmp = [data_all2randtmp; eeg_mi.trial{2}];
            
        end
        data_all1rand(:,:,rr) = data_all1randtmp((label_nr == roi_),:);
        data_all2rand(:,:,rr) = data_all2randtmp((label_nr == roi_),:);
        
    end
    %%
    n_ = numel(find((label_nr == roi_)));
    avg_1 = squeeze(nanmean(data_all1rand,1));
    avg_2 = squeeze(nanmean(data_all2rand,1));
    avg_D = squeeze(nanmean(data_all2rand-data_all1rand,1));
    
    
    pu1 = prctile(avg_1, 95, 2)';
    pl1 = prctile(avg_1, 5, 2)';
    av1 = nanmean(avg_1, 2)';
    
    pu2 = prctile(avg_2, 95, 2)';
    pl2 = prctile(avg_2, 5, 2)';
    av2 = nanmean(avg_2, 2)';
    
    
    puD = prctile((avg_2-avg_1), 95, 2)';
    plD = prctile((avg_2-avg_1), 5, 2)';
    avD = nanmean((avg_2-avg_1), 2)';
    
    run1avg = nanmean(data_all1(find(label_nr == roi_),:));
    run2avg =  nanmean(data_all2(find(label_nr == roi_),:));
    % get a p-value time course
    
    p_vals1 = sum(avg_1'>nanmean(data_all1(find(label_nr == roi_),:)))./1000;
    p_vals1(301) = 1;
    p_vals1(p_vals1==0) = 0.001;
    p_vals2 = sum(avg_2'>nanmean(data_all2(find(label_nr == roi_),:)))./1000;
    p_vals2(301) = 1;
    p_vals2(p_vals2==0) = 0.001;
    p_valsD = sum(avg_D'>nanmean(data_all2(find(label_nr == roi_),:)-...
        data_all1(find(label_nr == roi_),:)))./1000;
    p_valsD(301) = 1;
    p_valsD(p_valsD==0) = 0.001;
    
    tresh1 = fdr(p_vals1(150:451),0.05);
    tresh2 = fdr(p_vals2(150:451),0.05);
    treshD = fdr(p_valsD(150:451),0.05);
    
    
    if any(tresh1) || any(tresh2) || any(treshD) 
        figure;
        subplot(141);
        plot_ecog(ones(size(find(label_nr ==roi_))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos_all(find(label_nr ==roi_),:),...
            [], 0.2, [0 90], 22, 1, 2, repmat(yeo_colors(roi_+1,:), [256 1]));
        subplot(142);
        
        shadedErrorBar(eeg_mi.time{1}.*10, ...
            nanmean(data_all1(find(label_nr == roi_),:)) , ...
            nanstd(data_all1(find(label_nr == roi_),:))./sqrt(numel(find(label_nr == roi_))));
        hold on;
        
        shadedErrorBar(eeg_mi.time{1}.*10,...
             av1, [pu1-av1; -pl1+av1],'lineprops',{'color', [125 125 125]./255});
         
         xlabel('HC shift relative to peak (ms)');
         ylabel('Mutual Information');
         xlim([-1500 1500]);
%          ylim([0 yl_]);
         
         tmp = nan(size(eeg_mi.time{1}));
         tmp(p_vals1<tresh1) = 0 ;
         tmp(1:149)= nan;
         tmp(551:end)= nan;
         plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',1);
         
         title('run1')
         
         subplot(143);
         shadedErrorBar(eeg_mi.time{1}.*10, ...
             nanmean(data_all2(find(label_nr == roi_),:)) , ...
             nanstd(data_all2(find(label_nr == roi_),:))./sqrt(numel(find(label_nr == roi_))));
         
         hold on;
         shadedErrorBar(eeg_mi.time{1}.*10,...
             av2, [pu2-av2; -pl2+av2],'lineprops',{'color', [125 125 125]./255});
         
         xlabel('HC shift relative to peak (ms)');
         ylabel('Mutual Information');
         xlim([-1500 1500]);
         
         
         title('difference in HC vs. in Other (independent sample t-test')
%          ylim([0 yl_]);

         tmp = nan(size(eeg_mi.time{1}));
         tmp(p_vals2<tresh2) = 0 ;
         tmp(1:149)= nan;
         tmp(551:end)= nan;
         plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',1);
         
         title('run2');
         
         subplot(144);
         shadedErrorBar(eeg_mi.time{1}.*10, ...
             nanmean(data_all2(find(label_nr == roi_),:)- data_all1(find(label_nr == roi_),:)) , ...
             nanstd(data_all2(find(label_nr == roi_),:)-data_all1(find(label_nr == roi_),:))./sqrt(numel(find(label_nr == roi_))));
         
         
         hold on;
         shadedErrorBar(eeg_mi.time{1}.*10,...
             avD, [puD-avD; -plD+avD],'lineprops',{'color', [125 125 125]./255});
         
         
         xlabel('HC shift relative to peak (ms)');
         ylabel('Difference in Mutual Information');
         xlim([-1500 1500]);
         
         tmp = nan(size(eeg_mi.time{1}));
         tmp(p_valsD<treshD) = 0 ;
         tmp(1:149)= nan;
         tmp(551:end)= nan;
         plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',1);
         title('difference');
        
    end
    if roi_==18; 
        
        sig_all_shuff = (p_vals2<tresh2 & p_valsD<treshD) ;
        sig2_shuff = (p_vals2<tresh2) ;
    end
    
end
%% NOW do the big stats HC2>HC1 & HC2>OTHER2 & HC2-1>Other2-1 & HC2>HCshuffled & HC2-HC1 > HC2shuffled-HC1shuffled
roi_indices = label_nr == roi_; % roi is hippocampus, no other roi shows this effect!

% test HC vs. other again and fdr correct
[H, p2 , CI, stats] = ttest2(data_all2(find(roi_indices),:),...
    data_all2(find(~roi_indices&~eff_indices),:));
plim_2 = fdr(p2, 0.05);
[H, pdiff , CI, stats] = ttest2(data_all2(find(roi_indices),:) - ...
    data_all1(find(roi_indices),:), data_all2(find(~roi_indices&~eff_indices),:) - ...
    data_all1(find(~roi_indices&~eff_indices),:));
plim_d = fdr(pdiff, 0.05);

sig_all_other = p2< plim_2 & pdiff< plim_d;

[H, p , CI, stats] = ttest(data_all2(find(roi_indices),:),...
    data_all1(find(roi_indices),:));
thresh21 = fdr(p, 0.05);
sig21 = p<thresh21;


sig_strictest = sig21&sig_all_other&sig_all_shuff;

sig_strict2 = p2< plim_2 & sig2_shuff;
%% plot everything together with the corrected timepoints

figure;
subplot(241);
plot_ecog(ones(size(find(label_nr ==roi_))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos_all(find(label_nr ==roi_),:),...
    [], 0.2, [0 90], 22, 1, 2, repmat(yeo_colors(roi_+1,:), [256 1]));
subplot(245);
plot_ecog(ones(size(find(label_nr ==roi_))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos_all(find(label_nr ==roi_),:),...
    [], 0.2, [90 0], 22, 1, 2, repmat(yeo_colors(roi_+1,:), [256 1]));

subplot(2,4,[2;6]);

%plot shuffled
shadedErrorBar(eeg_mi.time{1}.*10,...
    av1, [pu1-av1; -pl1+av1],'lineprops',{'color', [207 218 230]./255});

hold on;

% plot other
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all1(find(~roi_indices&~eff_indices),:)) , ...
    nanstd(data_all1(find(~roi_indices&~eff_indices),:))./sqrt(numel(find(~roi_indices&~eff_indices))),...
    'lineprops',{'color', [40 156 174]./255});
hold on;
% plot hippocampus
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all1(find(label_nr == roi_),:)) , ...
    nanstd(data_all1(find(label_nr == roi_),:))./sqrt(numel(find(label_nr == roi_))),...
    'lineprops',{'color', [30 83 145]./255});

xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([-0.001 0.05]);

legend({'Shuffled' 'Other' 'Hippocampus'})
legend boxoff;
title('run 1');
set(gca,'fontname','calibri')


subplot(2,4,[3;7]);

% plot shuffled
shadedErrorBar(eeg_mi.time{1}.*10,...
    av2, [pu2-av2; -pl2+av2],'lineprops',{'color', [207 218 230]./255});

% plot other
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(~roi_indices&~eff_indices),:)) , ...
    nanstd(data_all2(find(~roi_indices&~eff_indices),:))./sqrt(numel(find(~roi_indices&~eff_indices))),...
    'lineprops',{'color', [222 172 139]./255});
hold on;

% plot hippocampus
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(label_nr == roi_),:)) , ...
    nanstd(data_all2(find(label_nr == roi_),:))./sqrt(numel(find(label_nr == roi_))),...
    'lineprops',{'color', [198 0 27]./255});

% plot significance
tmp = nan(size(eeg_mi.time{1}));
tmp(sig_strict2) = 0 ;
plot(eeg_mi.time{1}.*10, tmp, 'ok', 'MarkerFaceColor','k',  'MarkerSize',1);


xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-1500 1500]);
ylim([-0.001 0.05]);

legend({'Shuffled' 'Other' 'Hippocampus'})
legend boxoff;
title('run 2');
set(gca,'fontname','calibri')



subplot(2,4,[4;8]);

% plot shuffled
shadedErrorBar(eeg_mi.time{1}.*10,...
    avD, [puD-avD; -plD+avD],'lineprops',{'color',  [207 218 230]./255});


% plot other
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(~roi_indices&~eff_indices),:)- ...
    data_all1(find(~roi_indices&~eff_indices),:)) , ...
    nanstd(data_all2(find(~roi_indices&~eff_indices),:)-...
    data_all1(find(~roi_indices&~eff_indices),:))./...
    sqrt(numel(find(~roi_indices&~eff_indices)))...
    ,'lineprops',{'color',   [246 129 255]./255});

hold on;


% plot hippocampus
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(roi_indices),:)- ...
    data_all1(find(roi_indices),:)) , ...
    nanstd(data_all2(find(roi_indices),:)-...
    data_all1(find(roi_indices),:))./...
    sqrt(numel(find(roi_indices)))...
    ,'lineprops',{'color',  [116 78 144]./255});

% plot significance
tmp = nan(size(eeg_mi.time{1}));
tmp(sig_strictest) = -0.01;
plot(eeg_mi.time{1}.*10, tmp, 'ok', 'MarkerFaceColor','k',  'MarkerSize',1);


xlabel('HC shift relative to peak (ms)');
ylabel('Difference in Mutual Information');
xlim([-1500 1500]);
ylim([-0.011 0.025]);

legend({'Shuffled' 'Other' 'Hippocampus'})
legend boxoff;
title('Difference');
set(gca,'fontname','calibri')
%% which points are significant?
taxis = eeg_mi.time{1}.*10;

taxis(sig_strictest)
comfortable = taxis(p_valsD<0.05) 
%% NOW also plot what drives run2
prevec = sig_strict2;
postvec1 = sig_strict2;
postvec2 = sig_strict2;
prevec(301:end) = 0;
postvec1(1:301) = 0;
postvec1(345:end) = 0;
postvec2(1:345) = 0;

collectpre = zeros(4, 7); %raw delta theta alpha beta lowgamma highgamma
collectpost2 =  zeros(4, 7); %raw delta theta alpha beta lowgamma highgamma
collectpost2 =  zeros(4, 7); %raw delta theta alpha beta lowgamma highgamma

pre_all = nanmean(...
    nanmean(data_all2(find(label_nr == roi_),prevec)));

% get the average pre
collectpre(1,1) = nanmean(...
    nanmean(single_raw_all2(find(label_nr == roi_),prevec)));
collectpre(1,2) = nanmean(...
    nanmean(single_delta_all2(find(label_nr == roi_),prevec)));
collectpre(1,3) = nanmean(...
    nanmean(single_theta_all2(find(label_nr == roi_),prevec)));
collectpre(1,4) = nanmean(...
    nanmean(single_alpha_all2(find(label_nr == roi_),prevec)));
collectpre(1,5) = nanmean(...
    nanmean(single_beta_all2(find(label_nr == roi_),prevec)));
collectpre(1,6) = nanmean(...
    nanmean(single_gamma_low_all2(find(label_nr == roi_),prevec)));
collectpre(1,7) = nanmean(...
    nanmean(single_gamma_high_all2(find(label_nr == roi_),prevec)));

collectpre(2,1) = nanmean(...
    nanmean(drop_raw_all2(find(label_nr == roi_),prevec)))-pre_all;
collectpre(2,2) = nanmean(...
    nanmean(drop_delta_all2(find(label_nr == roi_),prevec)))-pre_all;
collectpre(2,3) = nanmean(...
    nanmean(drop_theta_all2(find(label_nr == roi_),prevec)))-pre_all;
collectpre(2,4) = nanmean(...
    nanmean(drop_alpha_all2(find(label_nr == roi_),prevec)))-pre_all;
collectpre(2,5) = nanmean(...
    nanmean(drop_beta_all2(find(label_nr == roi_),prevec)))-pre_all;
collectpre(2,6) = nanmean(...
    nanmean(drop_gamma_low_all2(find(label_nr == roi_),prevec)))-pre_all;
collectpre(2,7) = nanmean(...
    nanmean(drop_gamma_high_all2(find(label_nr == roi_),prevec)))-pre_all;

% get the std pre 

collectpre(3,1) = nanmean(...
    (1/n_).*nanstd(single_raw_all2(find(label_nr == roi_),prevec)));
collectpre(3,2) = nanmean(...
    (1/n_).*nanstd(single_delta_all2(find(label_nr == roi_),prevec)));
collectpre(3,3) = nanmean(...
    (1/n_).*nanstd(single_theta_all2(find(label_nr == roi_),prevec)));
collectpre(3,4) = nanmean(...
    (1/n_).*nanstd(single_alpha_all2(find(label_nr == roi_),prevec)));
collectpre(3,5) = nanmean(...
    (1/n_).*nanstd(single_beta_all2(find(label_nr == roi_),prevec)));
collectpre(3,6) = nanmean(...
    (1/n_).*nanstd(single_gamma_low_all2(find(label_nr == roi_),prevec)));
collectpre(3,7) = nanmean(...
    (1/n_).*nanstd(single_gamma_high_all2(find(label_nr == roi_),prevec)));

collectpre(4,1) = nanmean(...
    (1/n_).*nanstd(drop_raw_all2(find(label_nr == roi_),prevec)));
collectpre(4,2) = nanmean(...
    (1/n_).*nanstd(drop_delta_all2(find(label_nr == roi_),prevec)));
collectpre(4,3) = nanmean(...
    (1/n_).*nanstd(drop_theta_all2(find(label_nr == roi_),prevec)));
collectpre(4,4) = nanmean(...
    (1/n_).*nanstd(drop_alpha_all2(find(label_nr == roi_),prevec)));
collectpre(4,5) = nanmean(...
    (1/n_).*nanstd(drop_beta_all2(find(label_nr == roi_),prevec)));
collectpre(4,6) = nanmean(...
    (1/n_).*nanstd(drop_gamma_low_all2(find(label_nr == roi_),prevec)));
collectpre(4,7) = nanmean(...
    (1/n_).*nanstd(drop_gamma_high_all2(find(label_nr == roi_),prevec)));

% average first after ------------------------------------
collectpost1(1,1) = nanmean(...
    nanmean(single_raw_all2(find(label_nr == roi_),postvec1)));
collectpost1(1,2) = nanmean(...
    nanmean(single_delta_all2(find(label_nr == roi_),postvec1)));
collectpost1(1,3) = nanmean(...
    nanmean(single_theta_all2(find(label_nr == roi_),postvec1)));
collectpost1(1,4) = nanmean(...
    nanmean(single_alpha_all2(find(label_nr == roi_),postvec1)));
collectpost1(1,5) = nanmean(...
    nanmean(single_beta_all2(find(label_nr == roi_),postvec1)));
collectpost1(1,6) = nanmean(...
    nanmean(single_gamma_low_all2(find(label_nr == roi_),postvec1)));
collectpost1(1,7) = nanmean(...
    nanmean(single_gamma_high_all2(find(label_nr == roi_),postvec1)));

post1_all = nanmean(...
    nanmean(data_all2(find(label_nr == roi_),postvec1)));
collectpost1(2,1) = nanmean(...
    nanmean(drop_raw_all2(find(label_nr == roi_),postvec1)))-post1_all;
collectpost1(2,2) = nanmean(...
    nanmean(drop_delta_all2(find(label_nr == roi_),postvec1)))-post1_all;
collectpost1(2,3) = nanmean(...
    nanmean(drop_theta_all2(find(label_nr == roi_),postvec1)))-post1_all;
collectpost1(2,4) = nanmean(...
    nanmean(drop_alpha_all2(find(label_nr == roi_),postvec1)))-post1_all;
collectpost1(2,5) = nanmean(...
    nanmean(drop_beta_all2(find(label_nr == roi_),postvec1)))-post1_all;
collectpost1(2,6) = nanmean(...
    nanmean(drop_gamma_low_all2(find(label_nr == roi_),postvec1)))-post1_all;
collectpost1(2,7) = nanmean(...
    nanmean(drop_gamma_high_all2(find(label_nr == roi_),postvec1)))-post1_all;

% get first std after
collectpost1(3,1) = nanmean(...
    (1/n_).*nanstd(single_raw_all2(find(label_nr == roi_),postvec1)));
collectpost1(3,2) = nanmean(...
    (1/n_).*nanstd(single_delta_all2(find(label_nr == roi_),postvec1)));
collectpost1(3,3) = nanmean(...
    (1/n_).*nanstd(single_theta_all2(find(label_nr == roi_),postvec1)));
collectpost1(3,4) = nanmean(...
    (1/n_).*nanstd(single_alpha_all2(find(label_nr == roi_),postvec1)));
collectpost1(3,5) = nanmean(...
    (1/n_).*nanstd(single_beta_all2(find(label_nr == roi_),postvec1)));
collectpost1(3,6) = nanmean(...
    (1/n_).*nanstd(single_gamma_low_all2(find(label_nr == roi_),postvec1)));
collectpost1(3,7) = nanmean(...
    (1/n_).*nanstd(single_gamma_high_all2(find(label_nr == roi_),postvec1)));

collectpost1(4,1) = nanmean(...
    (1/n_).*nanstd(drop_raw_all2(find(label_nr == roi_),postvec1)));
collectpost1(4,2) = nanmean(...
    (1/n_).*nanstd(drop_delta_all2(find(label_nr == roi_),postvec1)));
collectpost1(4,3) = nanmean(...
    (1/n_).*nanstd(drop_theta_all2(find(label_nr == roi_),postvec1)));
collectpost1(4,4) = nanmean(...
    (1/n_).*nanstd(drop_alpha_all2(find(label_nr == roi_),postvec1)));
collectpost1(4,5) = nanmean(...
    (1/n_).*nanstd(drop_beta_all2(find(label_nr == roi_),postvec1)));
collectpost1(4,6) = nanmean(...
    (1/n_).*nanstd(drop_gamma_low_all2(find(label_nr == roi_),postvec1)));
collectpost1(4,7) = nanmean(...
    (1/n_).*nanstd(drop_gamma_high_all2(find(label_nr == roi_),postvec1)));


% average second after ------------------------------------

collectpost2(1,1) = nanmean(...
    nanmean(single_raw_all2(find(label_nr == roi_),postvec2)));
collectpost2(1,2) = nanmean(...
    nanmean(single_delta_all2(find(label_nr == roi_),postvec2)));
collectpost2(1,3) = nanmean(...
    nanmean(single_theta_all2(find(label_nr == roi_),postvec2)));
collectpost2(1,4) = nanmean(...
    nanmean(single_alpha_all2(find(label_nr == roi_),postvec2)));
collectpost2(1,5) = nanmean(...
    nanmean(single_beta_all2(find(label_nr == roi_),postvec2)));
collectpost2(1,6) = nanmean(...
    nanmean(single_gamma_low_all2(find(label_nr == roi_),postvec2)));
collectpost2(1,7) = nanmean(...
    nanmean(single_gamma_high_all2(find(label_nr == roi_),postvec2)));

post2_all = nanmean(...
    nanmean(data_all2(find(label_nr == roi_),postvec2)));
collectpost2(2,1) = nanmean(...
    nanmean(drop_raw_all2(find(label_nr == roi_),postvec2)))-post2_all;
collectpost2(2,2) = nanmean(...
    nanmean(drop_delta_all2(find(label_nr == roi_),postvec2)))-post2_all;
collectpost2(2,3) = nanmean(...
    nanmean(drop_theta_all2(find(label_nr == roi_),postvec2)))-post2_all;
collectpost2(2,4) = nanmean(...
    nanmean(drop_alpha_all2(find(label_nr == roi_),postvec2)))-post2_all;
collectpost2(2,5) = nanmean(...
    nanmean(drop_beta_all2(find(label_nr == roi_),postvec2)))-post2_all;
collectpost2(2,6) = nanmean(...
    nanmean(drop_gamma_low_all2(find(label_nr == roi_),postvec2)))-post2_all;
collectpost2(2,7) = nanmean(...
    nanmean(drop_gamma_high_all2(find(label_nr == roi_),postvec2)))-post2_all;

% std second after ------------------------------------

collectpost2(3,1) = nanmean(...
     (1/n_).*nanstd(single_raw_all2(find(label_nr == roi_),postvec2)));
collectpost2(3,2) = nanmean(...
     (1/n_).*nanstd(single_delta_all2(find(label_nr == roi_),postvec2)));
collectpost2(3,3) = nanmean(...
     (1/n_).*nanstd(single_theta_all2(find(label_nr == roi_),postvec2)));
collectpost2(3,4) = nanmean(...
     (1/n_).*nanstd(single_alpha_all2(find(label_nr == roi_),postvec2)));
collectpost2(3,5) = nanmean(...
     (1/n_).*nanstd(single_beta_all2(find(label_nr == roi_),postvec2)));
collectpost2(3,6) = nanmean(...
     (1/n_).*nanstd(single_gamma_low_all2(find(label_nr == roi_),postvec2)));
collectpost2(3,7) = nanmean(...
     (1/n_).*nanstd(single_gamma_high_all2(find(label_nr == roi_),postvec2)));

collectpost2(4,1) = nanmean(...
    (1/n_).*nanstd(drop_raw_all2(find(label_nr == roi_),postvec2)));
collectpost2(4,2) = nanmean(...
    (1/n_).*nanstd(drop_delta_all2(find(label_nr == roi_),postvec2)));
collectpost2(4,3) = nanmean(...
    (1/n_).*nanstd(drop_theta_all2(find(label_nr == roi_),postvec2)));
collectpost2(4,4) = nanmean(...
    (1/n_).*nanstd(drop_alpha_all2(find(label_nr == roi_),postvec2)));
collectpost2(4,5) = nanmean(...
    (1/n_).*nanstd(drop_beta_all2(find(label_nr == roi_),postvec2)));
collectpost2(4,6) = nanmean(...
    (1/n_).*nanstd(drop_gamma_low_all2(find(label_nr == roi_),postvec2)));
collectpost2(4,7) = nanmean(...
    (1/n_).*nanstd(drop_gamma_high_all2(find(label_nr == roi_),postvec2)));


%% plot it
bar(collectpre(1:2,:), 'stacked');
ylim([-0.05 0.05]);

figure
bar(collectpost1(1:2,:), 'stacked');
ylim([-0.05 0.05]);
figure

bar(collectpost2(1:2,:), 'stacked');
ylim([-0.05 0.05])
%%
subplot(121);
c = [0,0,0;distinguishable_colors(8, 'w')];
% 
% % plot hippocampus
% plot(eeg_mi.time{1}.*10, ...
%     nanmean(data_all2(find(label_nr == roi_),:)) -nanmean(nanmean(data_all2(find(label_nr == roi_),:))),...
%     'color',  c(1,:));
% hold on

% plot withour delta
l1 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_delta_all2(find(label_nr == roi_),:)) , ...
    'color', c(2,:));

hold on
% plot withour theta
l2 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_theta_all2(find(label_nr == roi_),:)) ,...
    'color', c(3,:));


hold on
% plot withour alpha
l3 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_alpha_all2(find(label_nr == roi_),:)) ,...
    'color', c(4,:));

hold on
% plot withour beta
l4 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_beta_all2(find(label_nr == roi_),:)) , ...
    'color', c(5,:));

hold on
% plot withour lg
l5 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_gamma_low_all2(find(label_nr == roi_),:)) ,...
    'color', c(6,:));

hold on
% plot withour hg
l6 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_gamma_high_all2(find(label_nr == roi_),:)) ,...
    'color', c(7,:));

hold on
% plot withour raw
l7 = plot(eeg_mi.time{1}.*10, ...
    nanmean(drop_raw_all2(find(label_nr == roi_),:)) ,...
    'color', c(8,:));



xlim([-1500 1500]);

% plot significance
tmp = nan(size(eeg_mi.time{1}));
tmp(sig_strictest) = 0.01;
plot(eeg_mi.time{1}.*10, tmp, 'ok', 'MarkerFaceColor','k',  'MarkerSize',3);
tmp = nan(size(eeg_mi.time{1}));
tmp(sig_strict2) = 0.011;
plot(eeg_mi.time{1}.*10, tmp, 'ob', 'MarkerFaceColor','b',  'MarkerSize',3);
ylim([0 0.04]);

hleg = legend([l1 l2 l3 l4 l5 l6 l7]);
hleg.String = ({'delta', 'theta', 'alpha', 'beta', 'lowGamma', 'highGamma', 'raw'});
legend boxoff
title('dropping each feature');

xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
%NOW also plot what drives run2 incrementing


% % plot hippocampus
% plot(eeg_mi.time{1}.*10, ...
%     nanmean(data_all2(find(label_nr == roi_),:)) -nanmean(nanmean(data_all2(find(label_nr == roi_),:))),...
%     'color',  c(1,:));
% hold on
subplot(122);
% plot withour delta
l1 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_delta_all2(find(label_nr == roi_),:)) , ...
    'color', c(2,:));

hold on
% plot withour theta
l2 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_theta_all2(find(label_nr == roi_),:)) ,...
    'color', c(3,:));


hold on
% plot withour alpha
l3 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_alpha_all2(find(label_nr == roi_),:)) ,...
    'color', c(4,:));

hold on
% plot withour beta
l4 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_beta_all2(find(label_nr == roi_),:)) , ...
    'color', c(5,:));

hold on
% plot withour lg
l5 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_gamma_low_all2(find(label_nr == roi_),:)) ,...
    'color', c(6,:));

hold on
% plot withour hg
l6 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_gamma_high_all2(find(label_nr == roi_),:)) ,...
    'color', c(7,:));

hold on
% plot withour raw
l7 = plot(eeg_mi.time{1}.*10, ...
    nanmean(single_raw_all2(find(label_nr == roi_),:)) ,...
    'color', c(8,:));



xlim([-1500 1500]);

% plot significance
tmp = nan(size(eeg_mi.time{1}));
tmp(sig_strictest) = -0.0003;
plot(eeg_mi.time{1}.*10, tmp, 'ok', 'MarkerFaceColor','k',  'MarkerSize',3);
tmp = nan(size(eeg_mi.time{1}));
tmp(sig_strict2) = 0;
plot(eeg_mi.time{1}.*10, tmp, 'ob', 'MarkerFaceColor','b',  'MarkerSize',3);
ylim([0 0.04]);
ylim([-0.002 0.012]);

hleg = legend([l1 l2 l3 l4 l5 l6 l7]);
hleg.String = ({'delta', 'theta', 'alpha', 'beta', 'lowGamma', 'highGamma', 'raw'});
legend boxoff
title('each feature alone');

xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');

%%




%% plot this roi
n_ = numel(find((label_nr == 1)));
yl_= 0.04;
% [H, p , CI, stats] = ttest(data_all2(find(hc_idices),:), ...
%     data_all1(find(hc_idices),:));
% plim_ = fdr(p, 0.05);

avg_1 = squeeze(nanmean(data_all1rand,1));
avg_2 = squeeze(nanmean(data_all2rand,1));

pu1 = prctile(avg_1, 95, 2)';
pl1 = prctile(avg_1, 5, 2)';
av1 = nanmean(avg_1, 2)';

pu2 = prctile(avg_2, 95, 2)';
pl2 = prctile(avg_2, 5, 2)';
av2 = nanmean(avg_2, 2)';


puD = prctile((avg_2-avg_1), 95, 2)';
plD = prctile((avg_2-avg_1), 5, 2)';
avD = nanmean((avg_2-avg_1), 2)';


figure;
subplot(131);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all1(find(hc_idices),:)) , ...
    nanstd(data_all1(find(hc_idices),:))./sqrt(n_));
hold on;

shadedErrorBar(eeg_mi.time{1}.*10,...
    av1, [pu1-av1; -pl1+av1],'lineprops',{'color', [125 125 125]./255});

xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-3000 3000]);
ylim([0 yl_]);

title('run1')

subplot(132);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(hc_idices),:)) , ...
    nanstd(data_all2(find(hc_idices),:))./sqrt(n_));

hold on;
shadedErrorBar(eeg_mi.time{1}.*10,...
    av2, [pu2-av2; -pl2+av2],'lineprops',{'color', [125 125 125]./255});

xlabel('HC shift relative to peak (ms)');
ylabel('Mutual Information');
xlim([-3000 3000]);
ylim([0 yl_]);
title('run2');

subplot(133);
shadedErrorBar(eeg_mi.time{1}.*10, ...
    nanmean(data_all2(find(hc_idices),:)- data_all1(find(hc_idices),:)) , ...
    nanstd(data_all2(find(hc_idices),:)-data_all1(find(hc_idices),:))./sqrt(n_));


hold on;
shadedErrorBar(eeg_mi.time{1}.*10,...
    avD, [puD-avD; -plD+avD],'lineprops',{'color', [125 125 125]./255});


xlabel('HC shift relative to peak (ms)');
ylabel('Difference in Mutual Information');
xlim([-3000 3000]);
ylim([- 0.02 0.02]);

%%

title('run2 - run1');

set(gcf, 'Position',  [600, 600, 1400, 400]);

hold on;
tmp = nan(size(eeg_mi.time{1}));
tmp(p<plim_) = 0 ;
plot(eeg_mi.time{1}.*10, tmp, 'or', 'MarkerFaceColor','r',  'MarkerSize',8);
