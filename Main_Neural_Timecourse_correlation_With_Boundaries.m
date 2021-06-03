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


%% select electrodes based on GMM fitting
electrode_selection_threshold = compute_electrode_selection_threshold;

%% post hoc, what drives the model uniquely? i.e. Build the model and project it
%thresh = 0.95;
all_proj_models1 = [];
all_proj_models2 = [];
all_As = {};
all_models = cell(size(sjs, 1),1);
all_dat = [];
all_epos = [];
n_chans = [];

for ss = 1 : size(sjs,1)
    ss;
    
    code_ = sjs(ss,:);
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));
    load(fullfile(deriv_pth, 'elec_orig.mat'));
    epos = getEposFromElec(elec_orig, model_eeg.label);
    
    all_models{ss} = model_eeg;
    all_dat = cat(1,all_dat, cell2mat(model_eeg.effect));
    all_epos = cat(1,all_epos, epos);
    n_chans = cat(1,n_chans, numel(model_eeg.label));
    
    model_eeg = all_models{ss};
    effect = (cell2mat(model_eeg.effect));
    
    select_vec = effect>= electrode_selection_threshold;
    
    As= model_eeg.A(select_vec);
    all_As{ss} = As;
    signal_of_1 = model_eeg.trial{1}(select_vec, :);
    signal_of_2 = model_eeg.trial{2}(select_vec, :);
    
    for aa = 1 : size(As,1)
        
        A = As{aa};
        
        % ---- compute 2 to 1 only 
        A2to1 = squeeze(A(1,2,:));
   
        p = length(A2to1); %modelorder
        p1 = p+1; %modelorder +1
        m = length(signal_of_2); 
        X = signal_of_2(aa,:)-mean( signal_of_2(aa,:)); %demeaned signal

        M = 1*(m-p); %n_variables *(length-lags) // vars is here 1
        %X0 = signal_of_1(aa, p1:m); not used
        
        %stack up the lags in a matrix
        XL = zeros(p,M);
        for k = 1:p
            XL(k,:) = X(:,p1-k:m-k,:); % concatenate trials for k-lagged observations  
        end
        
        model2to1_only = cat(2, zeros(1,p), A2to1'*XL);
        
        %project the model of 1 onto 1
        all_proj_models1 = [all_proj_models1; ...
            model2to1_only.*(signal_of_1(aa,:)-mean( signal_of_1(aa,:)))];
        % now compute 1 to 2 only        
        A1to2 = squeeze(A(2,1,:));
      
        p = length(A1to2); %modelorder
        p1 = p+1; %modelorder +1
        m = length(signal_of_1); 
        X = signal_of_1(aa,:)-mean( signal_of_1(aa,:)); %demeaned signal

        M = 1*(m-p); %n_variables *(length-lags) // vars is here 1
        %X0 = signal_of_1(aa, p1:m); not used
        
        %stack up the lags in a matrix
        XL = zeros(p,M);
        for k = 1:p
            XL(k,:) = X(:,p1-k:m-k,:); % concatenate trials for k-lagged observations  
        end
        
         
        model1to2_only = cat(2, zeros(1,p), A1to2'*XL);
        
           
        %project the model of 2 onto 2
        all_proj_models2 = [all_proj_models2; ...
            model1to2_only.*(signal_of_2(aa,:)-mean( signal_of_2(aa,:)))];
               
    end

end

%% smooth the projected models with average
all_proj_models1_sm = smoothdata(all_proj_models1,2,'movmean',100);
all_proj_models2_sm = smoothdata(all_proj_models2,2,'movmean',100);

%% Plot the projected time-course of prediction properly!
figure(1);
clf;
shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
    'lineprops',{'LineWidth',2, 'color', [131 131 255]./255,});
hold on;
set(gcf, 'color', 'w');
hold off;
%     ylabel('Brain -log(p)');
ylabel('Brain prediction (mean difference)');
xlabel('time (minutes)');
xlim([0 7.33])
set(gca,'fontsize', 12);
hold off;

set(gca,'fontsize', 18);
%% ======= START WITH THE BOUNDARIES =================================== %%

load([data_pth_ filesep 'sourcedata' filesep 'boundary_timecourses'])

proj_diff = all_proj_models1_sm - all_proj_models2_sm;
md = mean(proj_diff)';
%[H,P,CI,STATS] = ttest(all_proj_models1_sm, all_proj_models2_sm);
%
tmp1 = imresize(bound_tcs.tc1, [1 45000]);
tmp2 = imresize(bound_tcs.tc2, [1 45000]);
tc1 = tmp1(1:length(model_eeg.time{1}));
tc2 = tmp2(1:length(model_eeg.time{1}));
bth = (tmp1(1:length(model_eeg.time{1}))'+tmp2(1:length(model_eeg.time{1}))')./2;
dff = (tmp1(1:length(model_eeg.time{1}))'-tmp2(1:length(model_eeg.time{1}))');

%%
% the locking CPR to event boundaries

load([data_pth_ filesep 'sourcedata' filesep 'e2times205_sel.mat'])
bb_trials = [];
for bb = 1 : numel(e2times)
    bb_trials = [bb_trials;
        md(nearest(model_eeg.time{1}, e2times(bb)/1000)-600:...
        nearest(model_eeg.time{1}, e2times(bb)/1000)+600)'];
end
shadedErrorBar([-6000:10:6000], mean(bb_trials), std(bb_trials)./sqrt(19));

% now with random boundaries:
n_rand = 100000;
n_bounds =  19 %30 for appeal:
bb_trials_r = nan([size(bb_trials), n_rand]);
for rr = 1: n_rand
    rr
    e2timesr = randi([6010,10*length(md)-6010],n_bounds,1);
    for bb = 1 : numel(e2timesr)
        bb_trials_r(bb,:,rr) = md(nearest(model_eeg.time{1}, ...
            e2timesr(bb)/1000)-600:...
            nearest(model_eeg.time{1}, e2timesr(bb)/1000)+600)';
    end
end
rand_means = squeeze(mean(bb_trials_r))';
p_arr = (sum(abs(rand_means)> mean(bb_trials)))./n_rand;

shadedErrorBar([-6000:10:6000], mean(bb_trials), std(bb_trials)./sqrt(n_bounds));
hold on;
pu = prctile(squeeze(mean(bb_trials_r))', 95);
pl = prctile(squeeze(mean(bb_trials_r))', 5);
av = mean(squeeze(mean(bb_trials_r))');
shadedErrorBar([-6000:10:6000],...
     av, [pu-av; -pl+av],'lineprops',{'color', [207 218 230]./255});


%% cross correlation with avg. event boundaryness
n_rand = 1000;

[C_x, lags] = xcorr(md-mean(md),bth-mean(bth), 'coeff');

[mx_corr, mx_ind]  = max(C_x);
mx_lag_ms = lags(mx_ind)*10;

C_x_max_r = zeros(n_rand, 1);
C_xGrams_r = zeros(n_rand, length(C_x));
rand_md = squeeze(phaseran(md,n_rand));
for rr = 1 : n_rand
    rr
    mdr = rand_md(:,rr);
    [C_xr, lagsr] = xcorr(mdr-mean(mdr),bth-mean(bth), 'coeff');
    C_xGrams_r(rr,:) = C_xr;
    [C_x_max_r(rr), mx_indr]  = max(C_xr);
%     mx_lag_msr = lagsr(mx_indr)*10;
end
p_value = numel(find(C_x_max_r>max(C_x)))/n_rand
%% phase shuffling stats tc1
n_rand = 1000;

[C_x1, lags1] = xcorr(md-mean(md),tc1-mean(tc1), 'coeff');

[mx_corr, mx_ind]  = max(C_x1);
mx_lag_ms = lags1(mx_ind)*10;

C_x_max_r1 = zeros(n_rand, 1);
C_xGrams_r1 = zeros(n_rand, length(C_x1));
rand_md = squeeze(phaseran(md,n_rand));
for rr = 1 : n_rand
    rr
    mdr = rand_md(:,rr);
    [C_xr1, lagsr1] = xcorr(mdr-mean(mdr),tc1-mean(tc1), 'coeff');
    C_xGrams_r1(rr,:) = C_xr1;
    [C_x_max_r1(rr), mx_indr1]  = max(C_xr1);
%     mx_lag_msr = lagsr(mx_indr)*10;
end
%% phase shuffling stats tc2

[C_x2, lags2] = xcorr(md-mean(md),tc2-mean(tc2), 'coeff');

[mx_corr, mx_ind]  = max(C_x2);
mx_lag_ms = lags2(mx_ind)*10;

C_x_max_r2 = zeros(n_rand, 1);
C_xGrams_r2 = zeros(n_rand, length(C_x2));
rand_md = squeeze(phaseran(md,n_rand));
for rr = 1 : n_rand
    rr
    mdr = rand_md(:,rr);
    [C_xr2, lagsr2] = xcorr(mdr-mean(mdr),tc2-mean(tc2), 'coeff');
    C_xGrams_r2(rr,:) = C_xr2;
    [C_x_max_r2(rr), mx_indr2]  = max(C_xr1);
%     mx_lag_msr = lagsr(mx_indr)*10;
end
%% plot after phase shuffling

p_vals = sum(C_xGrams_r > C_x')./n_rand;
% fdr correct within +- 5 seconds
p_thresh = fdr(p_vals(floor(length(lags)/2 - 4000) : floor(length(lags)/2 + 4000)), 0.05);
sig_vec = p_vals < p_thresh;
ptc = nan(size(sig_vec));
ptc(sig_vec) = 0.275;

p_vals1 = sum(C_xGrams_r1 > C_x1')./n_rand;
% fdr correct within +- 5 seconds
p_thresh1 = fdr(p_vals1(floor(length(lags)/2 - 4000) : floor(length(lags)/2 + 4000)), 0.05);
sig_vec1 = p_vals1 < p_thresh1;
ptc1 = nan(size(sig_vec1));
ptc1(sig_vec1) = 0.265;


p_vals2 = sum(C_xGrams_r2 > C_x2')./n_rand;
% fdr correct within +- 5 seconds
p_thresh2 = fdr(p_vals2(floor(length(lags)/2 - 4000) : floor(length(lags)/2 + 4000)), 0.05);
sig_vec2 = p_vals2 < p_thresh2;
ptc2 = nan(size(sig_vec2));
ptc2(sig_vec2) = 0.255;


figure
plt_mean = mean(C_xGrams_r);
prz1u = prctile(C_xGrams_r,95,1) ;
prz1l = prctile(C_xGrams_r,5,1) ;

[C_x, lags] = xcorr(md-mean(md),bth-mean(bth), 'coeff');
p = plot(lags./100, C_x, 'LineWidth',2, 'color', [116 78 144]./255); 
hold on; 

shadedErrorBar(lagsr./100, plt_mean, [prz1u-plt_mean; -prz1l+plt_mean],'lineprops',{'color', [25 25 25]./255, 'linewidth', 2});


plt_mean = mean(C_xGrams_r1);
prz1u = prctile(C_xGrams_r1,95,1) ;
prz1l = prctile(C_xGrams_r1,5,1) ;

[C_x, lags] = xcorr(md-mean(md),tc1-mean(tc1), 'coeff');
p = plot(lags./100, C_x, 'LineWidth',2, 'color', [30 83 145]./255); 
hold on; 

shadedErrorBar(lagsr./100, plt_mean, [prz1u-plt_mean; -prz1l+plt_mean],'lineprops',{'color', [75 75 75]./255});

% p = line([0 0],[-0.25 0.3], 'color', 'black','LineStyle','--', 'LineWidth',2);
% set(p.Annotation.LegendInformation,'IconDisplayStyle','off');


plt_mean = mean(C_xGrams_r2);
prz1u = prctile(C_xGrams_r2,95,1) ;
prz1l = prctile(C_xGrams_r2,5,1) ;

[C_x, lags] = xcorr(md-mean(md),tc2-mean(tc2), 'coeff');
p = plot(lags./100, C_x, 'LineWidth',2, 'color', [198 0 27]./255); 
hold on; 

shadedErrorBar(lagsr./100, plt_mean, [prz1u-plt_mean; -prz1l+plt_mean],'lineprops',{'color', [125 125 125]./255});



p = plot(lags./100, ptc,  'o', 'MarkerEdgeColor', [116 78 144]./255, ...
    'MarkerFaceColor',[116 78 144]./255,  'MarkerSize',1);
 set(p.Annotation.LegendInformation,'IconDisplayStyle','off');

hold on; 

p1 = plot(lags./100, ptc1,  'o', 'MarkerEdgeColor', [30 83 145]./255, ...
    'MarkerFaceColor',[30 83 145]./255,  'MarkerSize',1);
 set(p1.Annotation.LegendInformation,'IconDisplayStyle','off');
 hold on; 

p2 = plot(lags./100, ptc2,  'o', 'MarkerEdgeColor', [198 0 27]./255, ...
    'MarkerFaceColor',[198 0 27]./255,  'MarkerSize',1);
 set(p2.Annotation.LegendInformation,'IconDisplayStyle','off');

 xlim([-4 4])
 ylim([-0.15 0.3]);
 legend({'both', 'shuffled','t1', 'shuffled', 't2', 'shuffled'})
 legend boxoff;
 ylabel('correlation');
 
 

%% Do statistics with permutation of run 1 and run 2 in the brain data
n_rand = 1000;
C_x_max_r = zeros(n_rand, 1);
for rr = 1 : n_rand
    rr
    rand_ =  ((rand(size(all_proj_models1,1),1) > 0.5));
    crand_ = ~rand_;
    modr1 = [diag(rand_),diag(crand_)] * [all_proj_models1_sm;all_proj_models2_sm];
    modr2 = [diag(crand_),diag(rand_)] * [all_proj_models1_sm;all_proj_models2_sm];
    
  %  [~,P,~,STATS] = ttest(modr1, modr2);
    mdr = mean(modr1-modr2)';
    [C_xtmp, lagstmp] = xcorr(mdr-mean(mdr),bth-mean(bth), 'coeff');
    [C_x_max_r(rr), ~] = max(C_xtmp);
end
[C_x, lags] = xcorr(md-mean(md),bth-mean(bth), 'coeff');
p_Value = numel(find(C_x_max_r>max(C_x)))/n_rand
%%   Plot the cross-correlogramm between brain prediction and event boundaryness (1 and 2 avg)
[C_x, lags] = xcorr(md-mean(md),bth-mean(bth), 'coeff');
plot(lags./100, C_x, 'LineWidth',2)

xlim([-75 75]);
xlabel('LAG (seconds)');
ylabel('correlation');

[C_x_max, Lag_maxI] = max(C_x); 
Lag_max = lags(Lag_maxI);
set(gcf, 'color', 'white');
set(gca,'fontsize', 18);
line([0 0],[-0.25 0.3], 'color', 'black','LineStyle','--', 'LineWidth',2);
ylim([-0.15 0.3])
%% plot the average boundaryness together with the prediction timecourse

figure('Renderer', 'painters', 'Position', [100 100 1200 400])
subplot(211)
plot(model_eeg.time{1}./60, bth, 'color', [0 114 189]./255, 'LineWidth',2);
hold on
ylabel('agreement avg run1/run2');
xlabel('minutes');
xlim([0 7.33])
set(gca,'fontsize', 12);
hold off

subplot(212);
shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
    'lineprops',{'color', [131 131 255]./255, 'LineWidth',2});
hold on

ylabel('Brain prediction (mean difference)');
xlabel('minutes');
xlim([0 7.33])
set(gca,'fontsize', 12);
hold off

set(gcf, 'color', 'w');

%% get p-value and plot histogram
p_val = numel(find(C_x_max_r > C_x_max))/n_rand
nhist(C_x_max_r,'color',  [0 114 189]./255);
set(gcf, 'color', 'w');
line([C_x_max C_x_max],[0 120], 'color', 'red','LineStyle','--', 'LineWidth',2);
xlabel('max correlation');
ylabel('number');
%xlim([-C_x_max-0.05 C_x_max+0.05])
set(gca,'fontsize', 12);

fprintf(['the maximal correlation is: ' num2str(C_x_max)]);
fprintf(['\n p =  ' num2str(p_val)]);
%% compare the first and second run of boundaryness in correlation with the brain-data
reject_vec = [23 50 51 88 89 91 92 93 94 104 112 156 165];
load([data_pth_ filesep 'sourcedata' filesep 'boundary_all_tc1']);
load([data_pth_ filesep 'sourcedata' filesep 'boundary_all_tc2']);
dts1(reject_vec,:) = [];
dts2(reject_vec, :) = [];
tc1 = zeros(size(dts1,1), 45000);
tc2 = zeros(size(dts2,1), 45000);
for ss = 1:size(dts1,1)
    ss
    tc1(ss,:) = imresize(dts1(ss,:), [1 45000]);
    tc2(ss,:)  = imresize(dts2(ss,:), [1 45000]);
end
tc1 = tc1(:, 1:length(model_eeg.time{1}));
tc2 = tc2(:, 1:length(model_eeg.time{1}));
%
d1 = mean(tc1);
d2 = mean(tc2);

[C_x1, lags] = xcorr(md-mean(md),d1'-mean(d1), 'coeff');
[C_x2, lags] = xcorr(md-mean(md),d2'-mean(d2), 'coeff');
real_diff = max(C_x2) - max(C_x1);

rand_N = 1000;
rand_diffs = zeros(rand_N, 1);
dts_both = cat(1, tc1, tc2);
for rr = 1 : rand_N
    rr
    dts_both = dts_both(randperm(size(dts_both, 1))',:);
    dts2r = dts_both(1:size(dts2,1), :);
    dts1r = dts_both(size(dts1,1)+1 : end, :);
    dt1 = mean(dts1r);
    dt2 = mean(dts2r);
    [C1, LAG] = xcorr(dt1'-mean(dt1), md-mean(md), 'coeff');
    [C2, LAG] = xcorr(dt2'-mean(dt2), md-mean(md), 'coeff');
    rand_diffs(rr) = max(C2)- max(C1);
   
end
%% get p-value and plot histogram
p_val = numel(find(rand_diffs > real_diff))/rand_N
nhist(rand_diffs,'color',  [0 114 189]./255);
set(gcf, 'color', 'w');
line([real_diff real_diff],[0 120], 'color', 'red','LineStyle','--', 'LineWidth',2);
xlabel('max correlation difference');
ylabel('number');
%xlim([-C_x_max-0.05 C_x_max+0.05])
set(gca,'fontsize', 12);
fprintf(['the maximal correlation difference is: ' num2str(real_diff)]);
fprintf(['\n p =  ' num2str(p_val)]);
%% CHECK: Do the cross correlation and phase-shuffeing electrode-wise and plot the max
n_rand = 1000;
randomize_ = true;
C_x_e = zeros(size(proj_diff,1),1);
C_x_p = zeros(size(proj_diff,1),1);
C_x_lags = zeros(size(proj_diff,1),1);
C_xGrams_e = zeros(size(proj_diff,1), length(C_x));
for ee = 1 : size(proj_diff,1)
    ee
    sig = proj_diff(ee,:);
    [C_x, lags] = xcorr(sig-mean(sig),bth-mean(bth), 'coeff');
    C_xGrams_e(ee,:) = C_x;
    [C_x_e(ee), m_idx] = max(C_x);
    C_x_lags(ee) = lags(m_idx)*10;
    
    
    C_x_max_r = zeros(n_rand, 1);
    if randomize_
        rand_e = squeeze(phaseran(sig',n_rand));
        for rr = 1 : n_rand
            %rr
            mdr = rand_e(:,rr);
            [C_xr, lagsr] = xcorr(mdr-mean(mdr),bth-mean(bth), 'coeff');
            [C_x_max_r(rr), mx_indr]  = max(C_xr);
            %     mx_lag_msr = lagsr(mx_indr)*10;
        end
        C_x_p(ee) = numel(find(C_x_max_r > C_x_e(ee)))/n_rand;
    end
end
%% this is to do the thing ROI based:fist check the correlation matrix for clusters...
tmp = triu(corr(proj_diff'))-eye(size(proj_diff,1));
tmp2 = tmp(tmp~=0);
[srt, inds] = sort(tmp(:,1));
tmp2 = corr(proj_diff(inds,:)');

%%

epos = all_epos(all_dat>=electrode_selection_threshold,:);
effect = all_dat(all_dat>=electrode_selection_threshold);

close all
figure;
plot_ecog(effect(C_x_p<0.05), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos(C_x_p<0.05,:),...
    [-0.003 0.003], 0.6, [-90 0], 22, 1, 2, []);
set(gcf, 'color', 'white')
%%
close all
figure;
plot_ecog(C_x_lags(C_x_p<0.05)./1000, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos(C_x_p<0.05,:),...
    [], 0.6, [-90 0], 22, 1, 2, jet(256));

view([-90 0]);
%%
%% now get the ROIS for the effect epos

epos = all_epos(all_dat>=electrode_selection_threshold,:);
label_nr = get_YEO_labels(epos);

% plot them in the YEO colors:

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

u_labels = unique(label_nr);
figure;
plot_ecog(ones(size(find(label_nr ==u_labels(1)))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos(find(label_nr ==u_labels(1)),:),...
    [], 0.2, [-75 15], 22, 1, 2, repmat(yeo_colors(u_labels(1)+1,:), [256 1]));
hold on
for ll = 2: length(u_labels)
    plot_ecog(ones(size(find(label_nr ==u_labels(ll)))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos(find(label_nr ==u_labels(ll)),:),...
    [], 0.2, [-75  15], 22, 0, 2, repmat(yeo_colors(u_labels(ll)+1,:), [256 1]));
end
set(gcf, 'color', 'white')

