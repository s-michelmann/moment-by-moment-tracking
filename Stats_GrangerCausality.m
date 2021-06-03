
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

%% PLOTTING The Brain in Subject Specific colors first:
alpha_ = 0.3;
colors = distinguishable_colors(size(sjs, 1), [1,1,1; 0,0,0;150/255,75/255,0] );
left_right = 2;
for ss = 1 : 9
    
    %subject code
    code_ = sjs(ss,:);
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    
    % load the original elec file
    load(fullfile(deriv_pth, 'elec_orig.mat'));
    
    %load the eeg to get the correct labels
    load([deriv_pth filesep 'eeg_manualica_notch.mat']);
    
    % get the correct electrode positions
    epos = getEposFromElec(elec_orig, eeg_manualica.label);
   
    if ss == 1     
        plot_ecog(ones(size(epos,1),1), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos,[],alpha_, [0 90], 22, 1, left_right, repmat(colors(ss,:), [256 1]));
    else
         plot_ecog(ones(size(epos,1),1), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos,[], alpha_, [0 90], 22, 0, left_right, repmat(colors(ss,:), [256 1]));
    end
end
set(gcf, 'color', 'white');
%%
cd('/Users/sm61/Google Drive/Projects/Pieman/Figures/electrode_placement')
%%
if left_right == 0
    print(gcf,'right_top.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'left_top.png','-dpng','-r300');
else
    print(gcf,'both_top.png','-dpng','-r300');
end
%%

view([-90 0]);
if left_right == 0
    print(gcf,'right_in.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'left_out.png','-dpng','-r300');
else
    print(gcf,'both_left.png','-dpng','-r300');
end
%%
view([90 0]);
if left_right ==0
    print(gcf,'right_out.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'left_in.png','-dpng','-r300');
else
    print(gcf,'both_right.png','-dpng','-r300');
end

%% plot the patient colors for the figure
figure; hold on;
tmp = nan(9,1);
for ss = 1 :9
    tmp(ss) = 1;
    plot([1:9], tmp, 'o', 'Markersize', 10, 'MarkerEdgeColor', colors(ss,:), 'MarkerFaceColor', colors(ss,:));
    hold on
    tmp(ss) = nan;
end
xlim([-1 10]);ylim([0 2]);
legend({'patient 1' ,'patient 2','patient 3','patient 4','patient 5',...
    'patient 6','patient 7','patient 8','patient 9'});
legend boxoff;
axis off;
set(gcf, 'color', 'w')

%%
close all
cd ([data_pth_ filesep 'code'])
%% start with the intuition of audio entraining the brain but not vice versa
view_ = [0 90];
left_right = 2;
alpha_ = 0.3;
vrtxsz = 22;
lim_ = [- 0.005 0.005]; % only for plotting together
plot_audio_to_brain = 0;
all_aud_dat = [];
for ss = 1 : size(sjs, 1)
    
    %subject code
    code_ = sjs(ss,:);
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];   
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];   

    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));   
    load(fullfile(deriv_pth, 'elec_orig.mat'));   
    [epos] = getEposFromElec(elec_orig, model_eeg.label);
    
    abs_ = zeros(numel(model_eeg.label), 1);
    bas_ = zeros(numel(model_eeg.label), 1);
    for ll = 1: numel(model_eeg.label)
        %run1 + run2  to brain from audio by 2
        abs_(ll) = (model_eeg.F{ll}(1,3) + model_eeg.F{ll}(2,3))/2;
        % run1 + run2  to audio from brain by 2
        bas_(ll) = (model_eeg.F{ll}(3,1) + model_eeg.F{ll}(3,2))/2;
    end
    all_aud_dat = cat(1, all_aud_dat, abs_-bas_);
    switch plot_audio_to_brain
        case 1
            if ss == 1;
                plot_ecog(abs_, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                    ,epos,lim_,alpha_, view_, vrtxsz, 1);
            else
                plot_ecog(abs_, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                    ,epos,lim_,alpha_, view_, vrtxsz, 0);
            end
        case -1
            if ss == 1;
                plot_ecog(bas_, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                    ,epos,lim_,alpha_, view_, vrtxsz, 1);
            else
                plot_ecog(bas_, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                    ,epos,lim_,alpha_, view_, vrtxsz, 0);
            end
        case 0
            
            if ss == 1;
                plot_ecog(abs_ - bas_, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                    ,epos,lim_,alpha_, view_, vrtxsz, 1, left_right);
            else
                plot_ecog(abs_ - bas_, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                    ,epos,lim_,alpha_, view_, vrtxsz, 0, left_right);
            end
    end
end

%%
cd('/Users/sm61/Google Drive/Projects/Pieman/Figures/audio_entrainment')
if left_right == 0
    print(gcf,'right_top.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'left_top.png','-dpng','-r300');
else
    print(gcf,'both_top.png','-dpng','-r300');
end
%%

view([-90 0]);
if left_right == 0
    print(gcf,'right_in.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'left_out.png','-dpng','-r300');
else
    print(gcf,'both_left.png','-dpng','-r300');
end
%%
view([90 0]);
if left_right ==0
    print(gcf,'right_out.png','-dpng','-r300');
elseif left_right == 1
    print(gcf,'left_in.png','-dpng','-r300');
else
    print(gcf,'both_right.png','-dpng','-r300');
end
%% save the colorbar as SVG
clf;
colorbar;
caxis([lim_]);
axis off

%%
cd ([data_pth_ filesep 'code']);

%% STATISTICS
% electrodewise first
n_rand = 1000;
dt = all_aud_dat;
rands_e = zeros(n_rand,1);
for rr = 1 : n_rand
    rand_ = [];
    rand_ =  ((rand(length(dt),1) > 0.5)*2 - 1); 
    rands_e(rr) = mean(dt.*rand_);
end
% hist(rands_, 100)
% line([mean(dt), mean(dt)], [mean(dt)+0.0001, n_rand/20] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
p_val_audio_electrodeWise = numel(find(rands_e>mean(dt)))/numel(rands_e);
fprintf([num2str(n_rand)...
    'permutations, p - value audio entrainment electrode permutation: '...
    num2str(p_val_audio_electrodeWise)]);
fprintf('\n');
% now subjectwise permutation as well
n_rand = 1000;
dt = all_aud_dat;
rands_s = zeros(n_rand,1);
for rr = 1 : n_rand
    rand_ = [];
    for cc = 1 : length(n_chans)
        rand_ = [rand_; ...
            ones(n_chans(cc),1)* ((rand(1,1) > 0.5)*2 - 1)];
    end
    rands_s(rr) = mean(dt.*rand_);
end
p_val_audio_subjectwise = numel(find(rands_s>mean(dt)))/numel(rands_s);
fprintf([num2str(n_rand)...
    'permutations, p - value audio entrainment subject permutation: '...
    num2str(p_val_audio_subjectwise)]);
fprintf('\n');
%% now pretty plot those results
close all
figure;
Violin(rands_e, 1, 'ViolinColor', [207 218 230 ]./255, 'Width', 0.3);
hold on;
plot(1.5, mean(dt), 'o', 'MarkerEdgeColor', [198 0 27]./255,...
    'MarkerFaceColor', [198 0 27]./255);
Violin(rands_s, 2, 'ViolinColor', [40 156 174]./255, 'Width', 0.3);
ylim([-0.0003 0.0003])
set(gca, 'xtick', [1:2]);
set(gca, 'xticklabels', {'electrode permutation', 'subject permutation'});
set(gcf, 'color', 'white');
set(gca, 'ytick', [-0.0003:0.0001:0.0003]);

set(gcf, 'Position', [100 100 400 1000]);
set(gca,'fontsize', 14);
set(gca,'fontname','calibri');
%% clear old stuff
clear dt rands_e rands_s all_aud_dat abs_ bas_ epos
%% =======================================================================%
%% Now the main analysis Granger Causality run2 minus run 1

%%
close all
cd ([data_pth_ filesep 'code']);
view_ = [0 90];
alpha_ = 0.3;
vrtxsz = 22;

threshold = false;
left_right_both = 2;

mesh_pth_ = '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/';
lim_ = [- 0.003 0.003]; % only for plotting together

lim_individual = false;

bar_on = false;
plot_together = true;
all_dat = [];
all_epos = [];
all_models = cell(size(sjs, 1),1);
n_chans = [];
sp = 1;
for ss = 1 : size(sjs, 1)
   % if ss == 8; continue; end %to test the stats without sj 8;
      %subject code
    code_ = sjs(ss,:);
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];    
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    %load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max.mat']));
    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));
    load(fullfile(deriv_pth, 'elec_orig.mat'));
    epos = getEposFromElec(elec_orig, model_eeg.label);
    
    all_dat = cat(1,all_dat, cell2mat(model_eeg.effect));
    all_epos = cat(1,all_epos, epos);
    n_chans = cat(1,n_chans, numel(model_eeg.label));
    all_models{ss} = model_eeg;
    
    if lim_individual
        lim_ = [-max(abs_((cell2mat(model_eeg.effect)))) ...
            max(abs_((cell2mat(model_eeg.effect))))];
    end
    if threshold
        if thresh
            select_vec = cell2mat(model_eeg.effect)>=thresh;
        end
        epos = epos(select_vec,:);
        plot_vec = plot_vec(select_vec,:);
    end

    if ~plot_together
        
        subplot(3,6,sp);
        
        %plot ecog left
        plot_ecog(cell2mat(model_eeg.effect), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos,lim_, alpha_, [90 0], vrtxsz, 1, 1,[], bar_on);
       % title(code_);
        sp = sp+1;
        %plot ecor right
         subplot(3,6,sp);
       
        %plot ecog left
        plot_ecog(cell2mat(model_eeg.effect), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos,lim_, alpha_, [-90 0], vrtxsz, 1, 0,[], bar_on);
        sp = sp+1;
    elseif plot_together  && ss == 1;
        plot_ecog(cell2mat(model_eeg.effect), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos,lim_, alpha_, view_, vrtxsz, 1, left_right_both,[], bar_on);
    else
        plot_ecog(cell2mat(model_eeg.effect), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
            ,epos,lim_, alpha_, view_ , vrtxsz, 0, left_right_both, [], bar_on);
    end  
end
%%
cd('/Users/sm61/Google Drive/Projects/Pieman/Figures/granger_causality_memory')
if left_right_both == 1
    print(gcf,'right_top.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_top.png','-dpng','-r300');
else
    print(gcf,'both_top.png','-dpng','-r300');
end
%%

view([-90 0]);
if left_right_both == 1
    print(gcf,'right_in.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_out.png','-dpng','-r300');
else
    print(gcf,'both_left.png','-dpng','-r300');
end
%%
view([90 0]);
if left_right_both ==1
    print(gcf,'right_out.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_in.png','-dpng','-r300');
else
    print(gcf,'both_right.png','-dpng','-r300');
end
%% save the colorbar as SVG
clf;
colorbar;
caxis([lim_]);
axis off
%% STATISTICS
% electrodewise first
n_rand = 1000;
dt = all_dat;
rands_e = zeros(n_rand,1);
for rr = 1 : n_rand
    rand_ = [];
    rand_ =  ((rand(length(dt),1) > 0.5)*2 - 1); 
    rands_e(rr) = mean(dt.*rand_);
end
% hist(rands_, 100)
% line([mean(dt), mean(dt)], [mean(dt)+0.0001, n_rand/20] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
p_val_memory_electrodeWise = numel(find(rands_e>mean(dt)))/numel(rands_e);
fprintf([num2str(n_rand)...
    'permutations, p - value memory electrode permutation: '...
    num2str(p_val_memory_electrodeWise)]);
fprintf('\n');
% now subjectwise permutation as well
n_rand = 1000;
dt = all_dat;
rands_s = zeros(n_rand,1);
for rr = 1 : n_rand
    rand_ = [];
    for cc = 1 : length(n_chans)
        rand_ = [rand_; ...
            ones(n_chans(cc),1)* ((rand(1,1) > 0.5)*2 - 1)];
    end
    rands_s(rr) = mean(dt.*rand_);
end
p_val_memory_subjectwise = numel(find(rands_s>mean(dt)))/numel(rands_s);
fprintf([num2str(n_rand)...
    'permutations, p - value memory subject permutation: '...
    num2str(p_val_memory_subjectwise)]);
fprintf('\n');
%% now pretty plot those results
close all
figure;
Violin(rands_e, 1, 'ViolinColor', [222 127 139 ]./255, 'Width', 0.3);
hold on;
plot(1.5, mean(dt), 'o', 'MarkerEdgeColor', [198 0 27]./255,...
    'MarkerFaceColor', [198 0 27]./255);
Violin(rands_s, 2, 'ViolinColor', [116 78 114]./255, 'Width', 0.3);
ylim([-0.0001 0.0001])
set(gca, 'xtick', [1:2]);
set(gca, 'xticklabels', {'electrode permutation', 'subject permutation'});
set(gcf, 'color', 'white');
set(gca, 'ytick', [-0.0001:0.0001:0.0001]);

set(gcf, 'Position', [100 100 400 1000]);
set(gca,'fontsize', 14);
ylabel('average difference (F)')
set(gca,'fontname','calibri');


%% now check if something is learned about the audio:
all_eff2 = [];
for ss = 1 : size(sjs,1)
    ss
    model_eeg = all_models{ss};
    Fs= model_eeg.F;
    effect2 = zeros(size(model_eeg.effect));
    for aa = 1 : size(Fs,1)
        
        
        F = Fs{aa};
        effect2(aa) = F(2,3) - F(1,3);        
    end
    all_eff2 = [all_eff2;effect2];
end

plot_ecog(all_eff2, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
                ,all_epos,[-0.0031 0.0031]);
%% STATISTICS
% electrodewise 
n_rand = 1000;
dt = all_eff2;
rands_ = zeros(n_rand,1);
for rr = 1 : n_rand
    rand_ = [];
    rand_ =  ((rand(length(dt),1) > 0.5)*2 - 1); 
    rands_(rr) = mean(dt.*rand_);
end
% hist(rands_, 100)
% line([mean(dt), mean(dt)], [mean(dt)+0.0001, n_rand/20] , 'LineStyle', '--', 'LineWidth', 2, 'Color', [1 0 1]);
p_val_audioLeaned_electrodeWise = numel(find(rands_<mean(dt)))/numel(rands_)
figure;
addpath('nhist');
nhist(rands_,'color', [173 174 249 ]./255);

%nhist(rands_,'color', [234 96 90]./255);
set(gcf, 'color', 'w');
line([mean(dt) mean(dt)],[0 2500], 'color', 'red','LineStyle','--', 'LineWidth',2);
xlabel('avererage prediction increase');
% xlim([min(rands_) mean(dt)+0.00001])
set(gca,'fontsize', 12);

%% now select electrodes based on GMM fitting
electrode_selection_threshold = compute_electrode_selection_threshold;
set(gca,'fontsize', 14);

set(gca,'fontname','calibri');
xlim([-0.002 0.004]);
%% now plot all electrodes in their corresponding yeo colors:

% these colors come from freesurfer
% plot ROIS in the YEO colors:
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

%% get the labels for each electrode position and the unique labels
label_nr = get_YEO_labels(all_epos);
u_labels = unique(label_nr);
%%
%% plot the yeo colors for the figure
close all
figure; hold on;
tmp = nan(size(yeo_colors,1),1);
for cc = 1 : size(yeo_colors,1)
    tmp(cc) = 1;
    plot([1:18], tmp, 'o', 'Markersize', 10, 'MarkerEdgeColor', ...
        yeo_colors(cc,:), 'MarkerFaceColor', yeo_colors(cc,:));
    hold on
    tmp(cc) = nan;
end
xlim([-1 19]);ylim([0 2]);
legend({'Yeo 0' ,'Yeo 1' ,'Yeo 2','Yeo 3','Yeo 4','Yeo 5',...
    'Yeo 6','Yeo 7','Yeo 8','Yeo 9','Yeo 10','Yeo 11' ,'Yeo 12','Yeo 13','Yeo 14','Yeo 15',...
    'Yeo 16','Yeo 17'});
legend boxoff;
axis off;
set(gcf, 'color', 'w')

%% plot the electrodes in yeo colors
cd([data_pth_ filesep 'code'])
close all;
view_ = [0 90];
alpha_ = 0.3;
vrtxsz = 22;
left_right_both = 2;

plot_ecog(ones(size(find(label_nr ==u_labels(1)))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,all_epos((label_nr ==u_labels(1)),:),...
    [], alpha_, view_, vrtxsz, 1, left_right_both, repmat(yeo_colors(u_labels(1)+1,:), [256 1]));
hold on
for ll = 2: length(u_labels)
    plot_ecog(ones(size(find(label_nr ==u_labels(ll)))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,all_epos((label_nr ==u_labels(ll)),:),...
    [], alpha_, view_, vrtxsz, 0, left_right_both, repmat(yeo_colors(u_labels(ll)+1,:), [256 1]));
end
set(gcf, 'color', 'white')
%%
cd('/Users/sm61/Google Drive/Projects/Pieman/Figures/yeo_labels_all')
if left_right_both == 1
    print(gcf,'right_top.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_top.png','-dpng','-r300');
else
    print(gcf,'both_top.png','-dpng','-r300');
end
%%

view([-90 0]);
if left_right_both == 1
    print(gcf,'right_in.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_out.png','-dpng','-r300');
else
    print(gcf,'both_left.png','-dpng','-r300');
end
%%
view([90 0]);
if left_right_both ==1
    print(gcf,'right_out.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_in.png','-dpng','-r300');
else
    print(gcf,'both_right.png','-dpng','-r300');
end
%% plot the effect electrodes in yeo colors
epos = all_epos(all_dat>=electrode_selection_threshold,:);
label_nr = get_YEO_labels(epos);

u_labels = unique(label_nr);
% count electrodes per yeo-roi
counts_ = arrayfun(@(x) [x numel(find(x == label_nr))],u_labels, 'Un', 0);
%%

cd([data_pth_ filesep 'code'])
close all;
view_ = [0 90];
alpha_ = 0.3;
vrtxsz = 22;
left_right_both = 1;
plot_ecog(ones(size(find(label_nr ==u_labels(1)))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos((label_nr ==u_labels(1)),:),...
    [], alpha_, view_, vrtxsz, 1, left_right_both, repmat(yeo_colors(u_labels(1)+1,:), [256 1]));
hold on
for ll = 2: length(u_labels)
    plot_ecog(ones(size(find(label_nr ==u_labels(ll)))), '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
    ,epos((label_nr ==u_labels(ll)),:),...
    [], alpha_, view_, vrtxsz, 0, left_right_both, repmat(yeo_colors(u_labels(ll)+1,:), [256 1]));
end
set(gcf, 'color', 'white')

%%
cd('/Users/sm61/Google Drive/Projects/Pieman/Figures/yeo_effect_labels')
if left_right_both == 1
    print(gcf,'right_top.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_top.png','-dpng','-r300');
else
    print(gcf,'both_top.png','-dpng','-r300');
end

view([-90 0]);
if left_right_both == 1
    print(gcf,'right_in.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_out.png','-dpng','-r300');
else
    print(gcf,'both_left.png','-dpng','-r300');
end
view([90 0]);
if left_right_both ==1
    print(gcf,'right_out.png','-dpng','-r300');
elseif left_right_both == 0
    print(gcf,'left_in.png','-dpng','-r300');
else
    print(gcf,'both_right.png','-dpng','-r300');
end
%%
cd([data_pth_ filesep 'code'])

%% post hoc, what drives the model uniquely? i.e. Build the model and project it
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
    %coefficient matrices
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

%% Plot the model order as a boxplot in YEO colors
m_orders = 10.*cell2mat(cellfun(@(x) size(x, 3), cat(1,all_As{:}), 'Un', 0));
%raincloud_plot(m_orders, [1 0 0 ]); hold on

bp = boxplot(m_orders,'positions', 1,'colors', [0 0 0], 'Widths',0.02);
set(bp,'linewidth',2);
h=findobj(gca,'tag','Outliers');
set(h,'Visible', 'off')
hold on
plot(1, mean(m_orders), 'o','Color', [198 0 27]./255, 'MarkerSize',7, 'MarkerFaceColor', [198 0 27]./255);

hold on;
for mm = 1:size(m_orders,1)
    plot(1.05+0.05*rand(1,1)-0.015, m_orders(mm), 'o','Color',...
        yeo_colors(label_nr(mm),:), 'MarkerSize',5, 'MarkerFaceColor', ...
        yeo_colors(label_nr(mm),:));
end
ylabel('model order (ms)');
set(gca,'xtick',[]);
xlim([0.9 1.15]);
ylim([100 400]);
set(gcf, 'color', 'white');
set(gca, 'ytick', [100:100:400]);

set(gcf, 'Position', [100 100 400 1000]);
set(gca,'fontsize', 18);
set(gca,'fontname','calibri');
%% also plot the lags as the absolute beta values
collect2to1 = nan(size(m_orders,1), max(m_orders)/10); 
collect1to2 = nan(size(m_orders,1), max(m_orders)/10); 
ii = 1;
for ss = 1 : numel(all_As)
    As =  all_As{ss};
    for aa = 1 : size(As,1)
        collect2to1(ii,end+1-size(As{aa},3):end) = (squeeze(As{aa}(1,2,:))');
        collect1to2(ii,end+1-size(As{aa},3):end) = (squeeze(As{aa}(2,1,:))');
        ii = ii +1;
    end
end
figure; 
plot([-350:10:-10], ((abs(collect2to1)-abs(collect1to2))),'Color',[0 0 0]+0.05*k);
hold on
plot([-350:10:-10], (nanmean((abs(collect2to1)-abs(collect1to2)))),'k')
xlabel('lag ms');
ylabel('beta value');
legend({'electrodes' 'mean'});
legend boxoff;
set(gca,'fontsize', 14);
set(gcf, 'color', 'white');
%% compare to the beta values of the auroregressive models
collect1to1 = nan(size(m_orders,1), max(m_orders)/10); 
collect2to2 = nan(size(m_orders,1), max(m_orders)/10); 
ii = 1;
for ss = 1 : numel(all_As)
    As =  all_As{ss};
    for aa = 1 : size(As,1)
        collect1to1(ii,end+1-size(As{aa},3):end) = (squeeze(As{aa}(1,1,:))');
        collect2to2(ii,end+1-size(As{aa},3):end) = (squeeze(As{aa}(2,2,:))');
        ii = ii +1;
    end
end
figure; 
plot([-350:10:-10], ((collect1to1 + collect2to2)./2),'Color',[0 0 0]+0.05*k);
hold on
plot([-350:10:-10], (nanmean((collect1to1 + collect2to2)./2)),'k')
xlabel('lag ms');
ylabel('beta value');

set(gca,'fontsize', 14);
set(gcf, 'color', 'white');
% 
% %% plot the model order on the brain
% close all
% figure;
% plot_ecog(m_orders, '/Users/sm61/Documents/tools/fieldtrip/template/anatomy/'...
%     ,all_epos(all_dat>=electrode_selection_threshold,:),...
%     [min(m_orders) max(m_orders)], 0.6, [-90 0], 22, 1, 2, plasma(256));
% 
% view([-90 0])
%% smooth the projected models with average
all_proj_models1_sm = smoothdata(all_proj_models1,2,'movmean',100);
all_proj_models2_sm = smoothdata(all_proj_models2,2,'movmean',100);
proj_diff = all_proj_models1_sm- all_proj_models2_sm;
%% check if the correlation between ROIS is lower than the avg. correlation 
%within rois
corr_within_rois = [];
corr_between_rois = [];
roi_tcs = zeros(length(u_labels), length(proj_diff));
for ll = 1: length(u_labels)
    diff_roi = proj_diff(label_nr==u_labels(ll),:);
    if size(diff_roi,1)>1
        c_mat = corr(diff_roi');
        tmp = triu(c_mat)-eye(size(c_mat));
        corr_within_rois = cat(1,corr_within_rois, tmp(tmp~=0));
        roi_tcs(ll,:) = mean(diff_roi,1);
    else
        roi_tcs(ll,:) = diff_roi;
    end
    
end
c_mat = corr(roi_tcs');
tmp = triu(c_mat)-eye(size(c_mat));
corr_between_rois = tmp(tmp~=0);

%%
Violin(corr_within_rois, 1, 'ViolinColor', [0.2 0.2 0.2]);
hold on
Violin(corr_between_rois, 2, 'ViolinColor', [0.2 0.2 0.2]);
set(gca, 'xtick', [1:2]);
set(gca, 'xticklabels', {'correlation within ROIs', 'correlation between ROIS'});
set(gcf, 'color', 'white');

[H, p, ci, stats] = ttest2(corr_within_rois, corr_between_rois);
disp(p);
title('Correlation of projected model-prediction difference')

%% Plot the projected time-course of prediction properly!
figure(1);
clf;
shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
    'lineprops',{'LineWidth',2, 'color', [39 40 56]./255});
% shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
%     'lineprops',{'LineWidth',2, 'color', [30 38 145]./255});
hold on;
set(gcf, 'color', 'w');
hold off;
%     ylabel('Brain -log(p)');
ylabel('mean projected model difference');
xlabel('time in story (minutes)');
xlim([0 7.33])
set(gca,'fontsize', 12);
hold off;

set(gca,'fontsize', 18);


set(gcf, 'Position', [100 100 2800 400]);
%% derive a phase-shuffled time-course of prediction
all_proj_models1_rand = [];
all_proj_models2_rand = [];

n_rand = 500;
for ss = 1 : size(sjs,1)
    ss;    
    code_ = sjs(ss,:)
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
 
    effect = (cell2mat(all_models{ss}.effect));
    select_vec = effect>= electrode_selection_threshold;
    %coefficient matrices
    
    signal_of_1 = all_models{ss}.trial{1}(select_vec, :);
    signal_of_2 = all_models{ss}.trial{2}(select_vec, :);
    As =   all_As{ss};
    for aa = 1 : size(As,1)
        aa
        A = As{aa};  
        % ---- compute 2 to 1 only 
        A2to1 = squeeze(A(1,2,:));
   
        p = length(A2to1); %modelorder
        p1 = p+1; %modelorder +1
        m = length(signal_of_2); 
        X = signal_of_2(aa,:)-mean( signal_of_2(aa,:)); %demeaned signal
        
        Xr = permute(phaseran(squeeze(X)', n_rand), [2,1,3]);
        M = 1*(m-p); %n_variables *(length-lags) // vars is here 1
        %X0 = signal_of_1(aa, p1:m); not used
        
        %stack up the lags in a matrix
        XLr = zeros(p,M, n_rand);
        for k = 1:p
            XLr(k,:,:) = Xr(:,p1-k:m-k,:); % concatenate trials for k-lagged observations  
        end
        
        model2to1_only_r = zeros(1, size(XLr,2)+p, n_rand);
         disp('rr 1');
        for rr = 1 : n_rand
            
            model2to1_only_r(:,:,rr) = cat(2, zeros(1,p), A2to1'*XLr(:,:,rr));
        end
        
        Xr2 = Xr;
     
        % now compute 1 to 2 only        
        A1to2 = squeeze(A(2,1,:));
      
        p = length(A1to2); %modelorder
        p1 = p+1; %modelorder +1
        m = length(signal_of_1); 
        X = signal_of_1(aa,:)-mean( signal_of_1(aa,:)); %demeaned signal
        Xr = permute(phaseran(squeeze(X)', n_rand), [2,1,3]);
        M = 1*(m-p); %n_variables *(length-lags) // vars is here 1
        %X0 = signal_of_1(aa, p1:m); not used
        
        %stack up the lags in a matrix
        XLr = zeros(p,M, n_rand);
        for k = 1:p
            XLr(k,:,:) = Xr(:,p1-k:m-k,:); % concatenate trials for k-lagged observations  
        end
        
        model1to2_only_r = zeros(1, size(XLr,2)+p, n_rand);
        disp('rr 2');
        for rr = 1 : n_rand
            
            model1to2_only_r(:,:,rr) = cat(2, zeros(1,p), A1to2'*XLr(:,:,rr));
        end
        
        Xr1 = Xr;
        %project the model of 1 onto 1
        all_proj_models1_rand = [all_proj_models1_rand; ...
            model2to1_only_r.*Xr1];
        %project the model of 2 onto 2
        all_proj_models2_rand = [all_proj_models2_rand; ...
            model1to2_only_r.*Xr2];
               
    end

end

all_proj_models1_sm_rand = smoothdata(all_proj_models1_rand,2,'movmean',100);
all_proj_models2_sm_rand = smoothdata(all_proj_models2_rand,2,'movmean',100);

% %% now shuffle 10000 times:
% n_rand = 10000;
% proj_diff_rand1 = cat(2,squeeze(...
%     mean(all_proj_models1_sm_rand - all_proj_models2_sm_rand)), ...
%     mean(all_proj_models1_sm - all_proj_models2)');
% proj_diff_rand = zeros(n_rand, size(proj_diff,2));
% for rr = 1 : n_rand
%     rr
%     
%     proj_diff_rand(rr,:) = sign(rand(size(proj_diff_rand1,2),1)-0.5)'*...
%         proj_diff_rand1';
% end
%%
proj_diff_rand = mean(all_proj_models1_sm_rand-all_proj_models2_sm_rand);
av = squeeze(mean((proj_diff_rand),3));
pu = squeeze(prctile((proj_diff_rand),95,3)); 
pl = squeeze(prctile((proj_diff_rand),5,3));
%% plot with 95th percentiles of random
figure(1);
clf;

% shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
%     'lineprops',{'LineWidth',2, 'color', [30 38 145]./255});

shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
    'lineprops',{'LineWidth',2, 'color', [39 40 56]./255});
hold on;
shadedErrorBar(model_eeg.time{1}./60,...
    av, [pu-av; -pl+av],'lineprops',{'color', [207 218 230]./255});
     
set(gcf, 'color', 'w');
hold off;
%     ylabel('Brain -log(p)');
ylabel('mean projected model difference');
xlabel('time in story (minutes)');
xlim([0 7.33])
set(gca,'fontsize', 12);
hold off;

set(gca,'fontsize', 18);


set(gcf, 'Position', [100 100 2800 400]);
%%
addpath(pwd);
%% RENDER THE VIDEO(s)
md = mean(all_proj_models1_sm-all_proj_models2_sm);

figure('Renderer', 'painters', 'Position', [100 100 2800 400])
v = VideoWriter('granger_predictionFinal.avi');
v.FrameRate = 60;
open(v);
for j=1000/v.FrameRate:1000/v.FrameRate:length(model_eeg.time{1})*10
    i = round(j);
    
    ix = nearest(model_eeg.time{1}, i/1000);
    
    figure(1);
    clf;
    shadedErrorBar(model_eeg.time{1}./60,mean(all_proj_models1_sm-all_proj_models2_sm), std(all_proj_models1_sm-all_proj_models2_sm)/sqrt(size(all_proj_models1,1)),...
        'lineprops',{'LineWidth',2, 'color',[39 40 56]./255,});
    hold on;
 
%     plot(model_eeg.time{1}(ix)/60, md(ix), 'ko','MarkerFaceColor',[50, 50, 50]./255)
    
    line([model_eeg.time{1}(ix)/60 model_eeg.time{1}(ix)/60],...
        [-0.04 0.1], 'color', [222 127 139]./255,'LineStyle','--', 'LineWidth',2);
    
    set(gcf, 'color', 'w');
    hold off;
    %     ylabel('Brain -log(p)');
    ylabel('mean projected model difference');
    xlabel('time in story (minutes)');
    xlim([0 7.33])
    set(gca,'fontsize', 12);
    
    set(gca,'fontname','calibri');
    hold off
    
    videoFrame = getframe(figure(1)); %step(videoFReader);
    writeVideo(v,videoFrame);
end

close(v);
