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

%% get the electrode-effect selection threshold
effect_thresh = compute_electrode_selection_threshold;
close all;
%% FOR EVERY SUBJECT
eff_indices = [];
epos_all = [];
hc_idices = [];
allmaps1 = [];
allmaps2 = [];
hc_labels = {};
for ss = 1 : size(sjs, 1)
    if ss == 1|| ss == 3 || ss == 8; continue; end
    % get sj code
    code_ = sjs(ss,:)
    % find path names
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    deriv_pth3 = [data_pth_ 'derivatives' filesep 'mi', filesep ...
        'sub-' code_ ];
    % load the EEG
    load([deriv_pth filesep 'eeg_manualica_notch.mat']);
    eeg = eeg_manualica; clear eeg_manualica;

    % load the GC results
    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));
    % keep the indices of the electrodess that show an effect above thresh
    eff_idx = (cell2mat(model_eeg.effect) >= effect_thresh);
    if ~any(eff_idx) 
        disp(['no effect for sj ' code_ ])
        continue;
        
    end    
    eff_indices = cat(1, eff_indices, eff_idx);
    
    load([deriv_pth, '/elec_orig.mat']);
    epos = getEposFromElec(elec_orig, eeg.label);
    epos_all = cat(1, epos_all, epos);
    % find 
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
    
    %
    try
    % load in the locked MI
    load(fullfile(deriv_pth3, ['LagAndShift2DMIaround_EB1s_Run1slidingWindow.mat']));
     load(fullfile(deriv_pth3, ['LagAndShift2DMIaround_EB1s_Run2slidingWindow.mat']));
    catch
        disp(['no data for sj ' code_ ])
        continue;
    end
    
    %reorganize the data and cut off the trailing nonsense...
    tmp1 = nan([size(eeg_mi1.trial{1},1) ,501, numel(eeg_mi1.trial)]);
    tmp2 = nan([size(eeg_mi1.trial{1},1) ,501, numel(eeg_mi1.trial)]);
    for tt = 1: numel(eeg_mi1.trial)
       % tt
        tmp1(:,1:501,tt) =  eeg_mi1.trial{tt}(:,1:501);
        tmp2(:,1:501,tt) =  eeg_mi2.trial{tt}(:,1:501);
    end
    % permute to channel * shift * lag
    allmaps1 = cat(1, allmaps1, permute(tmp1, [1, 3, 2]));
    allmaps2 = cat(1, allmaps2, permute(tmp2, [1, 3, 2]));
end
%% plot the hippocampus and effect channels in the corresponding colors
mkdir('/Users/sm61/Google Drive/Projects/Pieman/Figures/hc_and_effectElecs');
%%
addpath(pwd);
close all;
alpha_ = 0.15;
colors = [30/255,83/255,145/255;198/255,0/255,27/255] ;
left_right =2
plot_ecog(ones(numel(find(eff_indices)),1), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos_all(find(eff_indices),:),[],alpha_,...
    [0 90], 22, 1, left_right, repmat(colors(1,:), [256 1]));

plot_ecog(ones(numel(find(hc_idices)),1), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos_all(find(hc_idices),:),[],alpha_,...
    [0 90], 22, 1, left_right, repmat(colors(2,:), [256 1]));

set(gcf, 'color', 'white');
%
cd('/Users/sm61/Google Drive/Projects/Pieman/Figures/hc_and_effectElecs');
%
if left_right == 0
    print(gcf,'left_top.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'right_top.png','-dpng','-r300');
else
    print(gcf,'both_top.png','-dpng','-r300');
end
%

view([-90 0]);
if left_right == 0
    print(gcf,'left_out.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'right_in.png','-dpng','-r300');
else
    print(gcf,'both_left.png','-dpng','-r300');
end
%
view([90 0]);
if left_right ==0
    print(gcf,'left_in.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'right_out.png','-dpng','-r300');
else
    print(gcf,'both_right.png','-dpng','-r300');
end

%%   ============= STATISTICS!!! 

% covid19 hack ... work from home
cd ''
load('allmaps1.mat');
load('label_nr6sj.mat');
load('hc_idices6sj.mat');
load('rand_hc1.mat');
load('eff_idcs');
load('rand_hc1.mat');

%% stats compare to visual: 
%% RUN 1 

%% stats compare to visual: 

% cut the respective data within a generous but reasonable range,
% ignore lags outside of 1.5 seconds anyways, bc was changed in some
% control analyses scripts
hc_data_run1 = allmaps1(find((hc_idices)),100:350,151:351);
rand_hc_data_run1 = permute(rand_hc1(:,51:251,100:350), [3,2,1]);
other_data_run1 = allmaps1(find(~(hc_idices|eff_idcs)), 100:350,151:351);
vis_data_run1 = allmaps1(find(label_nr==1 &~(eff_idcs)), 100:350,151:351);

%get the real cluster
tmpmap = squeeze(nanmean(hc_data_run1)-nanmean(vis_data_run1));
tmpmap(isnan(tmpmap)) = 0;
[H, p , CI, stats] = ttest2(hc_data_run1, vis_data_run1, 'tail', 'both', 'alpha', 0.05);
H(isnan(H)) = 0;
clus = bwlabeln(squeeze(H));
max_sum = -Inf;
ccm = -1;
for cc = 1 : numel(unique(clus))
    tmp = sum(tmpmap(clus==cc));
    if tmp > max_sum
        max_sum = tmp;
        ccm = cc
    end
end
mxsum_real = max_sum;
msk_real = squeeze(clus==ccm);
% get the random clusters
% imagesc([-1000:10:1000],[-2000:10:500],  msk2.*squeeze(nanmean(hc_data_run1).*H)); axis xy
  all_dat4rand = cat(1, hc_data_run1, vis_data_run1);% 
  mxsums = zeros(1000,1)
for rr = 1 : 1000
    rr
    dat_rand = all_dat4rand(randperm(size(all_dat4rand,1)),:,:);
    hcr = dat_rand(1:size(hc_data_run1,1),:,:);
    visr = dat_rand(size(hc_data_run1,1)+1:end,:,:);
    tmpmap = squeeze(nanmean(hcr)-nanmean(visr));
    tmpmap(isnan(tmpmap)) = 0;
    [H, p , CI, stats] = ttest2(hcr, visr, 'tail', 'both', 'alpha', 0.05);
    H(isnan(H)) = 0;
    clus = bwlabeln(squeeze(H));
    max_sum = -Inf;
    for cc = 1 : numel(unique(clus))
        tmp = sum(tmpmap(clus==cc));
        if tmp > max_sum
            max_sum = tmp;
        end
    end
    mxsums(rr) = max_sum;
end

%%
permutation_vis_HC1000k_bthTail_05 = [];
permutation_vis_HC1000k_bthTail_05.mxsums_rand = mxsums;
permutation_vis_HC1000k_bthTail_05.pval = numel(find(mxsums>=mxsum_real))./1000;
permutation_vis_HC1000k_bthTail_05.mxsum_real = mxsum_real;
permutation_vis_HC1000k_bthTail_05.taxis_EB = [-2000:10:500];
permutation_vis_HC1000k_bthTail_05.lag_axis = [-1000:10:1000];
permutation_vis_HC1000k_bthTail_05.clusmap = msk_real;
save permutation1_vis_HC1000k_bthTail_05 permutation_vis_HC1000k_bthTail_05;

%%
% now  compare to other channels resampled: 
% imagesc([-1000:10:1000],[-2000:10:500],  msk2.*squeeze(nanmean(hc_data_run1).*H)); axis xy

% get the real cluster
tmpmap = squeeze(nanmean(hc_data_run1)-nanmean(other_data_run1));
tmpmap(isnan(tmpmap)) = 0;
[H, p , CI, stats] = ttest2(hc_data_run1, other_data_run1, 'tail', 'both', 'alpha', 0.05);
H(isnan(H)) = 0;
max_sum = -Inf;
ccm = -1;
for cc = 1 : numel(unique(clus))
    tmp = sum(tmpmap(clus==cc));
    if tmp > max_sum
        max_sum = tmp;
        ccm = cc
    end
end
mxsum_real = max_sum;
msk_real = squeeze(clus==ccm);

all_dat4rand = other_data_run1;%
mxsums = zeros(1000,1)
for rr = 1 : 1000
    rr
    dat_rand = all_dat4rand(randperm(size(all_dat4rand,1)),:,:);
    hcr = dat_rand(1:size(hc_data_run1,1),:,:);
    othr = dat_rand(size(hc_data_run1,1)+1:end,:,:);
    tmpmap = squeeze(nanmean(hcr)-nanmean(othr));
    tmpmap(isnan(tmpmap)) = 0;
    [H, p , CI, stats] = ttest2(hcr, othr, 'tail', 'both', 'alpha', 0.05);
    H(isnan(H)) = 0;
    clus = bwlabeln(squeeze(H));
    max_sum = -Inf;
    for cc = 1 : numel(unique(clus))
        tmp = sum(tmpmap(clus==cc));
        if tmp > max_sum
            max_sum = tmp;
        end
    end
    mxsums(rr) = max_sum;
end
%%
permutation_other_HC1000k_bthTail_05 = [];
permutation_other_HC1000k_bthTail_05.mxsums_rand = mxsums;
permutation_other_HC1000k_bthTail_05.pval = numel(find(mxsums>=mxsum_real))./1000;
permutation_other_HC1000k_bthTail_05.mxsum_real = mxsum_real;
permutation_other_HC1000k_bthTail_05.taxis_EB = [-2000:10:500];
permutation_other_HC1000k_bthTail_05.lag_axis = [-1000:10:1000];
permutation_other_HC1000k_bthTail_05.clusmap = msk_real;
save permutation1_other_HC1000k_bthTail_05 permutation_other_HC1000k_bthTail_05;

%% RUN 2
load('allmaps2.mat');
load('rand_hc2.mat');
load('rand_hc2.mat');

% select reasonable interval for statistics
hc_data_run2 = allmaps2(find((hc_idices)),100:350,151:351);
rand_hc_data_run2 = permute(rand_hc2(:,51:251,100:350), [3,2,1]);
other_data_run2 = allmaps2(find(~(hc_idices|eff_idcs)), 100:350,151:351);
vis_data_run2 = allmaps2(find(label_nr==1 &~(eff_idcs)), 100:350,151:351);

% this is the shuffled mask
msk2 = sum(rand_hc_data_run2>squeeze(nanmean(hc_data_run2)),3)<50;

% compute the real difference map
tmpmap = squeeze(nanmean(hc_data_run2)-nanmean(vis_data_run2));
tmpmap(isnan(tmpmap)) = 0;
[H, p , CI, stats] = ttest2(hc_data_run2, vis_data_run2, 'tail', 'both', 'alpha', 0.05);
H(isnan(H)) = 0;
clus = bwlabeln(squeeze(H));
max_sum = -Inf;
ccm = -1;
for cc = 1 : numel(unique(clus))
    tmp = sum(tmpmap(clus==cc));
    if tmp > max_sum
        max_sum = tmp;
        ccm = cc
    end
end
mxsum_real = max_sum;
msk_real = squeeze(clus==ccm);
% compute the random maps 
% imagesc([-1000:10:1000],[-2000:10:500],  msk2.*squeeze(nanmean(hc_data_run2).*H)); axis xy
all_dat4rand = cat(1, hc_data_run2, vis_data_run2);%
mxsums = zeros(1000,1)
for rr = 1 : 1000
    rr
    dat_rand = all_dat4rand(randperm(size(all_dat4rand,1)),:,:);
    hcr = dat_rand(1:size(hc_data_run2,1),:,:);
    visr = dat_rand(size(hc_data_run2,1)+1:end,:,:);
    tmpmap = squeeze(nanmean(hcr)-nanmean(visr));
    tmpmap(isnan(tmpmap)) = 0;
    [H, p , CI, stats] = ttest2(hcr, visr, 'tail', 'both', 'alpha', 0.05);
    H(isnan(H)) = 0;
    clus = bwlabeln(squeeze(H));
    max_sum = -Inf;
    for cc = 1 : numel(unique(clus))
        tmp = sum(tmpmap(clus==cc));
        if tmp > max_sum
            max_sum = tmp;
        end
    end
    mxsums(rr) = max_sum;
end
%%
permutation_vis_HC1000k_bthTail_05 = [];
permutation_vis_HC1000k_bthTail_05.mxsums_rand = mxsums;
permutation_vis_HC1000k_bthTail_05.pval = numel(find(mxsums>=mxsum_real))./1000;
permutation_vis_HC1000k_bthTail_05.mxsum_real = mxsum_real;
permutation_vis_HC1000k_bthTail_05.taxis_EB = [-2000:10:500];
permutation_vis_HC1000k_bthTail_05.lag_axis = [-1000:10:1000];
permutation_vis_HC1000k_bthTail_05.clusmap = msk_real;
save permutation2_vis_HC1000k_bthTail_05 permutation_vis_HC1000k_bthTail_05;

%% stats compare to other channels resampled: 
% imagesc([-1000:10:1000],[-2000:10:500],  msk2.*squeeze(nanmean(hc_data_run1).*H)); axis xy
tmpmap = squeeze(nanmean(hc_data_run2)-nanmean(other_data_run2));
tmpmap(isnan(tmpmap)) = 0;
[H, p , CI, stats] = ttest2(hc_data_run2, other_data_run2, 'tail', 'both', 'alpha', 0.05);
H(isnan(H)) = 0;
clus = bwlabeln(squeeze(H));
max_sum = -Inf;
ccm = -1;
for cc = 1 : numel(unique(clus))
    tmp = sum(tmpmap(clus==cc));
    if tmp > max_sum
        max_sum = tmp;
        ccm = cc
    end
end
mxsum_real = max_sum;
msk_real = squeeze(clus==ccm);
all_dat4rand = other_data_run2;%
mxsums = zeros(1000,1)
for rr = 1 : 1000
    rr
    dat_rand = all_dat4rand(randperm(size(all_dat4rand,1)),:,:);
    hcr = dat_rand(1:size(hc_data_run2,1),:,:);
    othr = dat_rand(size(hc_data_run2,1)+1:end,:,:);
    tmpmap = squeeze(nanmean(hcr)-nanmean(othr));
    tmpmap(isnan(tmpmap)) = 0;
    [H, p , CI, stats] = ttest2(hcr, othr, 'tail', 'both', 'alpha', 0.05);
    H(isnan(H)) = 0;
    clus = bwlabeln(squeeze(H));
    max_sum = -Inf;
    for cc = 1 : numel(unique(clus))
        tmp = sum(tmpmap(clus==cc));
        if tmp > max_sum
            max_sum = tmp;
        end
    end
    mxsums(rr) = max_sum;
end

%%
permutation_other_HC1000k_bthTail_05 = [];
permutation_other_HC1000k_bthTail_05.mxsums_rand = mxsums;
permutation_other_HC1000k_bthTail_05.pval = numel(find(mxsums>=mxsum_real))./1000;
permutation_other_HC1000k_bthTail_05.mxsum_real = mxsum_real;
permutation_other_HC1000k_bthTail_05.taxis_EB = [-2000:10:500];
permutation_other_HC1000k_bthTail_05.lag_axis = [-1000:10:1000];
permutation_other_HC1000k_bthTail_05.clusmap = msk_real;
save permutation2_other_HC1000k_bthTail_05 permutation_other_HC1000k_bthTail_05;


%% END of stats now some ====>

%=====> plotting

load('allmaps1.mat');
load('allmaps2.mat')
load('label_nr6sj.mat');
load('hc_idices6sj.mat');
load('rand_hc1.mat');
load('eff_idcs');
load('rand_hc1.mat');

%%
plot_run = 1;
%%
%to plot run 1
addpath('code\Colormaps')
if plot_run == 1

real_map = squeeze(nanmean(allmaps1(find(hc_idices),1:451,101:end-100)));
rand_map1 = squeeze(prctile(rand_hc1(:,:,1:451), 5, 1))';
rand_map2 = squeeze(prctile(rand_hc1(:,:,1:451), 95, 1))';
end
% to plot run 2
if plot_run == 2
real_map = squeeze(nanmean(allmaps2(find(hc_idices),1:451,101:end-100)));
rand_map1 = squeeze(prctile(rand_hc2(:,:,1:451), 5, 1))';
rand_map2 = squeeze(prctile(rand_hc2(:,:,1:451), 95, 1))';
end
%% plot with a shuffled version for orientation?
opengl software
addpath('D:\Google Drive\Projects\Pieman\PiemanECoG\code\')
lag_axis = [-1500:10:1500];
bound_axis = [-3000:10:1500];
real_map(:,151) = 0;
[mx, xmx] = max(max(real_map,[],2));
[~, ymx] = max(max(real_map,[],1));
surf(lag_axis, bound_axis, real_map, 'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud');
colormap(inferno(256));
freezeColors()
hold on;
surf(lag_axis, bound_axis, rand_map1, 'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud', 'FaceAlpha', 0.6);
colormap(repmat([207 218 230]./255, 256, 1));
hold on;
surf(lag_axis, bound_axis, rand_map2, 'FaceColor','interp',...
   'EdgeColor','none',...
   'FaceLighting','gouraud', 'FaceAlpha', 0.6);
colormap(repmat([207 218 230]./255, 256, 1));

view([6 50]);
zlim([0 0.045] );
ylim([-3000 1000]);
set(gcf, 'color', 'white');
xlabel('lag');
%ylabel('time around boundary')
zlabel('MI');
plot3([lag_axis(ymx); lag_axis(ymx)],[bound_axis(xmx) bound_axis(xmx)],...
    [0 mx-0.01], 'k--');
hold on;
plot3([lag_axis(ymx); 1500],[bound_axis(xmx) bound_axis(xmx)],...
    [0 0], 'k--');
hold on;
plot3([lag_axis(ymx); lag_axis(ymx)],[-3000 bound_axis(xmx)],...
    [0 0], 'k--')
%%
load('permutation1_vis_HC1000k_bthTail_05_may27.mat');
%
figure;
Violin(permutation_vis_HC1000k_bthTail_05.mxsums_rand, 1, 'ViolinColor', [30 83 145 ]./255, 'Width', 0.3);
hold on;
plot(1.2, mean(permutation_vis_HC1000k_bthTail_05.mxsum_real), 'o', 'MarkerEdgeColor', [198 0 27]./255,...
    'MarkerFaceColor', [255 0 0]./255);


load('permutation2_vis_HC1000k_bthTail_05_may27.mat');
Violin(permutation_vis_HC1000k_bthTail_05.mxsums_rand, 2, 'ViolinColor',  [198 0 27]./255, 'Width', 0.3);
hold on;
plot(2.2, mean(permutation_vis_HC1000k_bthTail_05.mxsum_real), 'o', 'MarkerEdgeColor', [198 0 27]./255,...
    'MarkerFaceColor', [255 0 0]./255);


load('permutation1_other_HC1000k_bthTail_05_march30.mat');
Violin(permutation_other_HC1000k_bthTail_05.mxsums_rand, 3, 'ViolinColor', [40 156 174]./255, 'Width', 0.3);
hold on;
plot(3.2, mean(permutation_other_HC1000k_bthTail_05.mxsum_real), 'o', 'MarkerEdgeColor', [198 0 27]./255,...
    'MarkerFaceColor', [255 0 0]./255);
hold on;

load('permutation2_other_HC1000k_bthTail_05_march30.mat');
%
Violin(permutation_other_HC1000k_bthTail_05.mxsums_rand, 4, 'ViolinColor', [219 95 91]./255, 'Width', 0.3);
hold on;
plot(4.2, mean(permutation_other_HC1000k_bthTail_05.mxsum_real), 'o', 'MarkerEdgeColor', [198 0 27]./255,...
    'MarkerFaceColor', [255 0 0]./255);


set(gca, 'xtick', [1:4]);
set(gca, 'xticklabels', {'visual control', '', 'random other',''});
set(gcf, 'color', 'white');
%set(gca, 'ytick', [-0.0003:0.0001:0.0003]);

set(gcf, 'Position', [100 100 200 600]);
set(gca,'fontsize', 14);
set(gca,'fontname','calibri');
%% unmasked visual and HC for the supplemental information .. change this should be the main figure!


real_map1 = squeeze(nanmean(allmaps1(find(hc_idices),1:451,101:end-100)));

vis_map1 = squeeze(nanmean(allmaps1(find(label_nr==1),1:451,101:end-100)));

real_map2 = squeeze(nanmean(allmaps2(find(hc_idices),1:451,101:end-100)));

vis_map2 = squeeze(nanmean(allmaps2(find(label_nr==1),1:451,101:end-100)));

%%
real_map = vis_map2;

lag_axis = [-1500:10:1500];
bound_axis = [-3000:10:1500];
real_map(:,151) = 0;
[mx, xmx] = max(max(real_map,[],2));
[~, ymx] = max(max(real_map,[],1));
imagesc(lag_axis, bound_axis, real_map); axis xy;
colormap(jet(256));
caxis([0 0.045]);

ylim([-3000 1000]);
set(gcf, 'color', 'white');
xlabel('lag');
ylabel('time around EB');
set(gcf, 'color', 'white');
colorbar

