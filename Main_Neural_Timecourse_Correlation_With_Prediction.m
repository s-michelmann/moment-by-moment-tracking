
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

%% build the model again and project it
%thresh = 0.95;
all_proj_models1 = [];
all_proj_models2 = [];
all_As = {};
all_models = cell(size(sjs, 1),1);
all_dat = [];
all_epos = [];
n_chans = [];
all_signals_1 =[];
all_signals_2 =[];

for ss = 1 : size(sjs,1)
    ss;
    
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
    
    all_signals_1 = cat(1, all_signals_1, signal_of_1);
    
    all_signals_2 = cat(1, all_signals_2, signal_of_2);
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

all_signals_1_sm  = smoothdata(all_signals_1,2,'movmean',100);
all_signals_2_sm  = smoothdata(all_signals_2,2,'movmean',100);

%% get the ROIS for the effect epos and plot them again in yeo colors
epos = all_epos(all_dat>=electrode_selection_threshold,:);
label_nr = get_YEO_labels(epos);

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

%% Relate the memory guided prediction to the cloze probability Amy's data
word_path_a = [data_pth_ filesep 'sourcedata' filesep  ...
    'word_prediction' filesep 'amy']; 
word_path_b = [data_pth_ filesep 'sourcedata' filesep  ...
    'word_prediction' filesep 'replication'];

proj_diff = all_proj_models1_sm - all_proj_models2_sm;
md = mean(proj_diff)';

load([word_path_a filesep 'cloze1_a']);
load([word_path_a filesep 'cloze2_a']);
load([word_path_a filesep 'words_a.mat']);

load([word_path_a filesep 'word_onsets_a.mat']);
load([word_path_a filesep 'word_offsets_a.mat']);

[h,p, ci, stats] = ttest(run2_cloze, run1_cloze)

% collect time course for cloze 1, cloze 2 and the difference
tc_cloze_1a = zeros(size(md));
tc_cloze_2a = zeros(size(md));
cloze_delta_a = run2_cloze-run1_cloze;

ind_list_a = zeros(size(onsets,1),length(model_eeg.time{1}));
for ww = 1 : size(onsets,1)
    ind_list_a(ww,:) = (model_eeg.time{1}>onsets(ww) & model_eeg.time{1}<offsets(ww));
    tc_cloze_1a(ind_list_a(ww,:)==1) = run1_cloze(ww);
    tc_cloze_2a(ind_list_a(ww,:)==1) = run2_cloze(ww);
end
tc_cloze_delta_a = tc_cloze_2a-tc_cloze_1a;

%% Relate the memory guided prediction to the cloze probability Bobbi's data

load([word_path_b filesep 'cloze1_b_exclude']);
load([word_path_b filesep 'cloze2_b_exclude']);
load([word_path_b filesep 'words_b.mat']);

load([word_path_b filesep 'word_onsets_b.mat']);
load([word_path_b filesep 'word_offsets_b.mat']);

% collect time course for cloze 1, cloze 2 and the difference

ind_list_b = zeros(size(onsets_b,1),length(model_eeg.time{1}));
tc_cloze_1b = zeros(size(md));
tc_cloze_2b = zeros(size(md));
cloze_delta_b = cloze2_b-cloze1_b;
for ww = 3 : size(onsets_b,1)
    ind_list_b(ww,:) = (model_eeg.time{1}>onsets_b(ww) & model_eeg.time{1}<offsets_b(ww));
    tc_cloze_1b(ind_list_b(ww,:)==1) = cloze1_b(ww);
    tc_cloze_2b(ind_list_b(ww,:)==1) = cloze2_b(ww);
end
tc_cloze_delta_b = tc_cloze_2b-tc_cloze_1b;

corr(tc_cloze_delta_b(tc_cloze_delta_a~=0&tc_cloze_delta_b~=0)...
    ,tc_cloze_delta_a(tc_cloze_delta_a~=0&tc_cloze_delta_b~=0))

%% Relate the memory guided prediction to the cloze probability Bobbi's data cosine similarity

load([word_path_b filesep 'individual_cosSimexclude.mat']);
load([word_path_b filesep 'individual_cosSim2exclude.mat']);
cosine_delta_b = nanmean(individual_cosSim2,2) - nanmean(individual_cosSim,2);
nan_inds = (isnan(cosine_delta_b));
words_b_nonan = words_b(~nan_inds);
onsets_b_nonan = onsets_b(~nan_inds);
offsets_b_nonan = offsets_b(~nan_inds);
individual_cosSim_nonan = individual_cosSim(~nan_inds,:);
individual_cosSim2_nonan = individual_cosSim2(~nan_inds,:);

cosine_delta_b_nonan = nanmean(individual_cosSim2_nonan,2) - nanmean(individual_cosSim_nonan,2);
ind_list_b_nonan = zeros(size(onsets_b_nonan,1),length(model_eeg.time{1}));

tc_cosine_1b = zeros(size(md));
tc_cosine_2b = zeros(size(md));
for ww = 1 : size(onsets_b_nonan,1)
    ind_list_b_nonan(ww,:) = (model_eeg.time{1}>onsets_b_nonan(ww) & model_eeg.time{1}<offsets_b_nonan(ww));   
    tc_cosine_1b(ind_list_b_nonan(ww,:)==1) = ...
        nanmean(individual_cosSim_nonan(ww,:),2);
    tc_cosine_2b(ind_list_b_nonan(ww,:)==1) = ...
        nanmean(individual_cosSim2_nonan(ww,:),2);
end
tc_cosine_delta_b = tc_cosine_2b-tc_cosine_1b;
%%
sec_lag = 2;
fs = 100;
slsamp = sec_lag*fs;
%% for amy's data
n_rand = 1000;

true_ca = zeros(size(-slsamp:1:slsamp));
rand_ca =  zeros(n_rand, length(-slsamp:1:slsamp));
sig = tc_cloze_delta_a;
% compute the true lagged correlation between neural prediction and cloze
% probability
for ll = -slsamp:1:slsamp
    ll
    tic;
    md_shift = circshift(md, -ll);
    true_ca(ll+slsamp+1) = corr(md_shift(sig~=0), ...
        sig(sig~=0));
end
% and collect random
for rr = 1: n_rand
    rr
    sig_rand = zeros(size(sig));
    rand_cloze_delta = cloze_delta_a(randperm(length(words)));
    for ww = 1 : size(onsets,1)
        sig_rand(ind_list_a(ww,:)==1) = rand_cloze_delta(ww);
    end
    for ll = -slsamp:1:slsamp
        
        md_shift = circshift(md, -ll);

        rand_ca(rr,ll+slsamp+1) = corr(md_shift(sig_rand~=0), ...
            sig_rand(sig_rand~=0));
    end  
end

%% for bobbi's data
true_cb = zeros(size(-slsamp:1:slsamp));
rand_cb =  zeros(n_rand, length(-slsamp:1:slsamp));
sig = tc_cloze_delta_b;
% compute the true lagged correlation between neural prediction and cloze
% probability
for ll = -slsamp:1:slsamp
    ll
    tic;
    md_shift = circshift(md, -ll);
    true_cb(ll+slsamp+1) = corr(md_shift(sig~=0), ...
        sig(sig~=0));
end
% and collect random
for rr = 1: n_rand
    rr
    sig_rand = zeros(size(sig));
    tmp = cloze_delta_b(randperm(length(words_b)));
    rand_cloze_delta = [nan; nan; tmp(~isnan(tmp))];
    for ww = 3 : size(onsets_b,1)
        sig_rand(ind_list_b(ww,:)==1) = rand_cloze_delta(ww);
    end
    for ll = -slsamp:1:slsamp        
        md_shift = circshift(md, -ll);
        rand_cb(rr,ll+slsamp+1) = corr(md_shift(sig_rand~=0), ...
            sig_rand(sig_rand~=0));
    end  
end
%% for bobbi's data - change in glove distance

true_cbcs = zeros(size(-slsamp:1:slsamp));
rand_cbcs =  zeros(n_rand, length(-slsamp:1:slsamp));
sig = tc_cosine_delta_b;

% compute the true lagged correlation between neural prediction and cloze
% probability
for ll = -slsamp:1:slsamp
    ll
    tic;
    md_shift = circshift(md, -ll);
    true_cbcs(ll+slsamp+1) = corr(md_shift(sig~=0), ...
        sig(sig~=0));
end
% and collect random
for rr = 1: n_rand
    rr
    sig_rand = zeros(size(sig));
    rand_cloze_delta = cosine_delta_b_nonan(randperm(length(words_b_nonan)));
    for ww = 1 : size(onsets_b_nonan,1)
        sig_rand(ind_list_b_nonan(ww,:)==1) = rand_cloze_delta(ww);
    end
    for ll = -slsamp:1:slsamp        
        md_shift = circshift(md, -ll);
        rand_cbcs(rr,ll+slsamp+1) = corr(md_shift(sig_rand~=0), ...
            sig_rand(sig_rand~=0));
    end  
end
%%
ava = squeeze(mean((rand_ca)));
pua = squeeze(prctile((rand_ca),95)); 
pla = squeeze(prctile((rand_ca),5));

avb = squeeze(mean((rand_cb)));
pub = squeeze(prctile((rand_cb),95)); 
plb = squeeze(prctile((rand_cb),5));

avbcs = squeeze(mean((rand_cbcs)));
pubcs = squeeze(prctile((rand_cbcs),95)); 
plbcs = squeeze(prctile((rand_cbcs),5));

p_amy = sum(rand_ca>true_ca)./n_rand;
sig_amy = p_amy<fdr(p_amy, .05);
p_bobbi = sum(rand_cb>true_cb)./n_rand;
sig_bobbi = p_bobbi<fdr(p_bobbi, .05);

p_bobbics = sum(rand_cbcs>true_cbcs)./n_rand;
sig_bobbics = p_bobbics<fdr(p_bobbics, .05);

%% make a pretty figure
figure; 

% plot amy's data
plot([-slsamp:1:slsamp].*1000/fs, true_ca, 'Linewidth', 2, 'color', [116 78 144]./255);
hold on;


% plot bobbi's data
plot([-slsamp:1:slsamp].*1000/fs, true_cb, 'Linewidth', 2, 'color',[198 0 27]./255);
hold on;


% plot bobbi's data2
plot([-slsamp:1:slsamp].*1000/fs, true_cbcs, 'Linewidth', 2, 'color', [219 95 91]./255);
hold on;

%%

shadedErrorBar( [-slsamp:1:slsamp].*1000/fs,...
    ava, [pua-ava; -pla+ava],'lineprops',{'Linewidth', 2,'color', [40 156 174  ]./255});
hold on;

shadedErrorBar( [-slsamp:1:slsamp].*1000/fs,...
    avb, [pub-avb; -plb+avb],'lineprops',{'Linewidth', 2,'color', [246 129 255]./255});
hold on;


ylim([-0.1 0.2]);
shadedErrorBar( [-slsamp:1:slsamp].*1000/fs,...
    avbcs, [pubcs-avbcs; -plbcs+avbcs],'lineprops',{'Linewidth', 2,'color', [39 40 56]./255});
hold on;


%make a vector with the significance points
plcv = nan(size(true_ca));
plcv(sig_amy) = 0.19;
% plot that vector
plot([-slsamp:1:slsamp].*1000/fs, plcv, 'o', 'MarkerEdgeColor', [116 78 144]./255, ...
    'MarkerFaceColor',[116 78 144]./255,  'MarkerSize',1);

%make a vector with the significance points
plcv = nan(size(true_cbcs));
plcv(sig_bobbi) = 0.185;

plot([-slsamp:1:slsamp].*1000/fs, plcv,  'o', 'MarkerEdgeColor', [198 0 27]./255, ...
    'MarkerFaceColor',[198 0 27]./255,  'MarkerSize',1);

%make a vector with the significance points
plcv = nan(size(true_cb));
plcv(sig_bobbics) = 0.18;

plot([-slsamp:1:slsamp].*1000/fs, plcv,  'o', 'MarkerEdgeColor', [219 95 91]./255, ...
    'MarkerFaceColor',[198 0 27]./255,  'MarkerSize',1);

legend({'exp 1 prob', 'exp 2 prob', 'exp 2 vec', 'shuffled 1',  'shuffled 2',  'shuffled 3' });
legend boxoff;
xlabel('time around word onset');
ylabel('correlation');
set(gcf, 'color', 'white');

set(gca,'fontname','calibri');

set(gca,'Fontsize',12);

set(gcf, 'Position',  [600, 600, 600, 800]);
