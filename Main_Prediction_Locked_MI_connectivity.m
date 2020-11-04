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
%% Loading Prediction RESULTS

all_dat = {};
all_models = cell(size(sjs, 1),1);
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
    all_dat = cat(1,all_dat, model_eeg.effect);
    all_models{ss} = model_eeg;
    
end
%% derive the model for selected electrodes

% this will collect the models from run2 predicting run1, projected on run1
sj_modelsProjon1 = cell(size(sjs,1),1);

% this will collect the models from run1 predicting run2, projected on run2
sj_modelsProjon2 = cell(size(sjs,1),1);

for ss = 1 : size(sjs,1)
    
    model_eeg = all_models{ss};
    effect = (cell2mat(model_eeg.effect));
    
    % which electrodes show an effect?
    select_vec = effect>= electrode_selection_threshold;
    
    % get the model coefficients for these electrodes
    As= model_eeg.A(select_vec);
    
    % get the gamma power at these electrodes
    signal_of_1 = model_eeg.trial{1}(select_vec, :);
    signal_of_2 = model_eeg.trial{2}(select_vec, :);
    
    % loop through the electrodes
    for aa = 1 : size(As,1)
        
        %pull out the coefficient matrix
        A = As{aa};
        
        % ---- compute 2 to 1 only
        A2to1 = squeeze(A(1,2,:));
        p = length(A2to1); %modelorder
        p1 = p+1; %modelorder +1
        m = length(signal_of_2);
        X = signal_of_2(aa,:)-mean( signal_of_2(aa,:)); %demeaned signal
        
        M = 1*(m-p); %n_variables *(length-lags) // vars is here 1
        
        %stack up the lags in a matrix
        XL = zeros(p,M);
        for k = 1:p
            XL(k,:) = X(:,p1-k:m-k,:); % concatenate trials for k-lagged observations
        end
        
        % predicted data is coefficients times stacked Matrix
        model2to1_only = cat(2, zeros(1,p), A2to1'*XL);
        
        %project the models of 1 onto the signals of 1
        sj_modelsProjon1{ss,1} =  [sj_modelsProjon1{ss,1} ; ...
            model2to1_only.*(signal_of_1(aa,:)-mean( signal_of_1(aa,:)))];
        
        % now compute 1 to 2 only
        A1to2 = squeeze(A(2,1,:));
        p = length(A1to2); %modelorder
        p1 = p+1; %modelorder +1
        m = length(signal_of_1);
        X = signal_of_1(aa,:)-mean( signal_of_1(aa,:)); %demeaned signal
        
        M = 1*(m-p); %n_variables *(length-lags) // vars is here 1
        
        %stack up the lags in a matrix
        XL = zeros(p,M);
        for k = 1:p
            XL(k,:) = X(:,p1-k:m-k,:); % concatenate trials for k-lagged observations
        end
        
        % predicted data is coefficients times stacked Matrix
        model1to2_only = cat(2, zeros(1,p), A1to2'*XL);
        
        %project the models of 2 onto the signals 2
        sj_modelsProjon2{ss,1} =  [sj_modelsProjon2{ss,1} ; ...
            model1to2_only.*(signal_of_2(aa,:)-mean( signal_of_2(aa,:)))];
    end
end

%% now compute the mutual information at peak-locked data

sampling_rate = all_models{1}.fsample;
%loop though subjects
for ss = 1 : size(sjs, 1)
    
    
    disp(['sj = ' num2str(ss)])
       
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
    load([deriv_pth filesep 'eeg_manualica_notch.mat']);
    eeg = eeg_manualica; clear eeg_manualica;
    
    % load in the hippocampal channels
    fclose all;
    fid = fopen([deriv_pth filesep 'HC_channels.tsv'], 'r');
    txt = textscan(fid, '%s', 'delimiter', '\t', 'Headerlines', 0);
    fclose all   ;
    
    % if there are no channels continue with the next subject
    if isempty(txt{1})
        continue;
    end
    
    % unnecessary delete later
    % keep the indices of the HC channels for later plotting?
    hc_idx = find( cell2mat(cellfun(@(x) ismember(x, txt{:}), eeg.label, 'Un', 0)));
    
    % load in the Granger Causality effects
    load(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']));
    
    % keep the indices of the electrodess that show an effect above thresh
    eff_idx = find(cell2mat(model_eeg.effect) >= electrode_selection_threshold);
    
    % also continue if there is no effect
    if isempty(eff_idx)
        continue;
    end
    
    %names of the channels that show an effect
    eff_chans = model_eeg.label(eff_idx);
    % select an eeg for the effect and for the hippocampus
    cfg = []; cfg.channel = eff_chans;
    eeg_effect = ft_selectdata(cfg, eeg);
    
    %% USE ALL CHANNELs
    cfg = []; cfg.channel = model_eeg.label;
    eeg_all = ft_selectdata(cfg, eeg);
    
    clear eeg;
    
    % bpfilter eeg in the high gamma band (35-55Hz)
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [70 200];
    eeg_hgall = ft_preprocessing(cfg, eeg_all);
    eeg_hgeff = ft_preprocessing(cfg, eeg_effect);
    
    % bpfilter eeg in the low gamma band (35-55Hz)
    cfg.bpfreq = [35 55];
    eeg_lgall = ft_preprocessing(cfg, eeg_all);
    
    % bpfilter eeg in the beta band (15-30Hz)
    cfg.bpfreq = [15 30];
    eeg_ball = ft_preprocessing(cfg, eeg_all);
    
    % bpfilter eeg in the alpha band (8-15Hz)
    cfg.bpfreq = [8 15];
    eeg_aall = ft_preprocessing(cfg, eeg_all);
    
    % bpfilter eeg in the theta band (4-8Hz)
    cfg.bpfreq = [4 8];
    eeg_tall = ft_preprocessing(cfg, eeg_all);
    
    % lpfilter eeg for the delta band (min-4Hz)
    cfg.bpfilter = 'no';
    cfg.lpfilter = 'yes';
    cfg.lpfreq = 4;
    eeg_dall = ft_preprocessing(cfg, eeg_all);
    
    
    % loop through all the channels to get power via absolute of the
    % hilbert transform                                                         (TODO; try this with phase as well!?)
    for cc = 1 : numel(eeg_all.label)
        % and replace the data with the amplitude of the hilbert
        eeg_hgall.trial{1}(cc, :) = abs(hilbert( eeg_hgall.trial{1}(cc, :)));
        eeg_hgall.trial{2}(cc, :) = abs(hilbert( eeg_hgall.trial{2}(cc, :)));
        
        eeg_lgall.trial{1}(cc, :) = abs(hilbert( eeg_lgall.trial{1}(cc, :)));
        eeg_lgall.trial{2}(cc, :) = abs(hilbert( eeg_lgall.trial{2}(cc, :)));
        
        eeg_ball.trial{1}(cc, :) = abs(hilbert( eeg_ball.trial{1}(cc, :)));
        eeg_ball.trial{2}(cc, :) = abs(hilbert( eeg_ball.trial{2}(cc, :)));
        
        eeg_aall.trial{1}(cc, :) = abs(hilbert( eeg_aall.trial{1}(cc, :)));
        eeg_aall.trial{2}(cc, :) = abs(hilbert( eeg_aall.trial{2}(cc, :)));
        
        eeg_tall.trial{1}(cc, :) = abs(hilbert( eeg_tall.trial{1}(cc, :)));
        eeg_tall.trial{2}(cc, :) = abs(hilbert( eeg_tall.trial{2}(cc, :)));
        
        eeg_dall.trial{1}(cc, :) = abs(hilbert( eeg_dall.trial{1}(cc, :)));
        eeg_dall.trial{2}(cc, :) = abs(hilbert( eeg_dall.trial{2}(cc, :)));
    end
    % loop through all the Effect channels to get power via absolute of the
    for cc = 1 : numel(eeg_effect.label)
        % and replace the data with the amplitude of the hilbert
        eeg_hgeff.trial{1}(cc, :) = abs(hilbert( eeg_hgeff.trial{1}(cc, :)));
        eeg_hgeff.trial{2}(cc, :) = abs(hilbert( eeg_hgeff.trial{2}(cc, :)));
    end
    
    % downsample to 100Hz
    cfg = [];
    cfg.resamplefs = 100;
    
    eeg_all = ft_resampledata(cfg, eeg_all);
    eeg_hgall = ft_resampledata(cfg, eeg_hgall);
    eeg_lgall = ft_resampledata(cfg, eeg_lgall);
    eeg_ball = ft_resampledata(cfg, eeg_ball);
    eeg_aall = ft_resampledata(cfg, eeg_aall);
    eeg_tall = ft_resampledata(cfg, eeg_tall);
    eeg_dall = ft_resampledata(cfg, eeg_dall);
    
    eeg_hgeff = ft_resampledata(cfg, eeg_hgeff);
    % finally detrend the data to remove drift (but mostly mean!)
    cfg = [];
    cfg.detrend = 'yes';
    eeg_all  = ft_preprocessing(cfg, eeg_all);
    eeg_hgall = ft_preprocessing(cfg, eeg_hgall);
    eeg_lgall = ft_preprocessing(cfg, eeg_lgall);
    eeg_ball = ft_preprocessing(cfg, eeg_ball);
    eeg_aall = ft_preprocessing(cfg, eeg_aall);
    eeg_tall = ft_preprocessing(cfg, eeg_tall);
    eeg_dall = ft_preprocessing(cfg, eeg_dall);
    
    eeg_hgeff = ft_preprocessing(cfg, eeg_hgeff);
    
    
    % define the window of possible lags (-3s to +3s)
    window_width = 6*eeg_hgeff.fsample;
    
    % define the window with to compute mutual information across
    window_width_small = eeg_hgeff.fsample;
    
    % create a dummy for the mutual information data
    eeg_mi = rmfield(eeg_hgall, 'cfg');
    eeg_mi.sampleinfo = [1 window_width+1; window_width+2 2*window_width+2];
    
    % for each trial initialize the trial with zero and set the time to the
    % sampling points
    for tt = 1 : 2
        eeg_mi.trial{tt} = zeros(numel(eeg_all.label), ...
            window_width+1);
        eeg_mi.time{tt} = [-window_width/2:100/eeg_all.fsample:window_width/2];
    end
   
    % the labels are all labels
    eeg_mi.label = eeg_all.label;
    %% now get the peaks
    %first smooth the projected models with a moving average of 1 second
    proj_modelsOn1_sm = smoothdata(sj_modelsProjon1{ss,1},2,'movmean',sampling_rate);
    proj_modelsOn2_sm = smoothdata(sj_modelsProjon2{ss,1},2,'movmean',sampling_rate);
    
    % get the difference between model for 1 and model for 2
    proj_diff = mean(proj_modelsOn1_sm - proj_modelsOn2_sm);
    
    
    
    %% get peaks and troughs
    
    % make sure that no peak is detected in the first 3 and last 3 seconds
    pad_sec = 3;
    
    %--- this is the peak detection procedure
    %smooth the projected difference with a 2s gaussian for peak detection
    %which highlights peaks
    smd = smoothdata(proj_diff,'gaussian',200);
    
    % set the first and last pad_seconds to 0
    smd(1:(sampling_rate*pad_sec)) = 0;
    smd(end-(sampling_rate*pad_sec-1):end) = 0;
    
    % only moments exceeding 95th percentile are interesting. Label those
    % intervals in blocks
    peak1 = bwlabeln(smd > prctile(smd, 95));
    
    % everything else is set to nan
    smd(peak1==0) = nan;
    
    %initialize a vector of zeros for peak-times
    ptimes = zeros(max(peak1),1);
    
    %loop through all blocks of high values and find the maxima
    for pp = 1: max(peak1)
        zvec = zeros(size(smd)); %initialize a zero_vector
        
        %set the zero vector to the actual values only in the block
        zvec(peak1 == pp) = smd(peak1 == pp);
        
        % find the maimum inside this block
        [~, ptimes(pp)] = max(zvec);
    end
    
    % get the time-points in the first run that correspond to these
    % sampling points
    ptimes_s = model_eeg.time{1}(ptimes);
    
    
    
    %%
    
    % loop through trials
    for trial_idx  =1 : 2
        
        % select this trial!
        cfg = [];
        cfg.trials = trial_idx;       
        eeg_hg2 = ft_selectdata(cfg, eeg_hgall);
        eeg_lg2 = ft_selectdata(cfg, eeg_lgall);
        eeg_b2 = ft_selectdata(cfg, eeg_ball);
        eeg_a2 = ft_selectdata(cfg, eeg_aall);
        eeg_t2 = ft_selectdata(cfg, eeg_tall);
        eeg_d2 = ft_selectdata(cfg, eeg_dall);
        eeg_2  = ft_selectdata(cfg, eeg_all);
        eeg_hgeff2 = ft_selectdata(cfg, eeg_hgeff);
        
        % loop through channels
        for cc = 1 : numel(eeg_hg2.label)
            
            fprintf('\nsj = %d trial = %d channel = %d',ss, trial_idx, cc);
            
            %initialize the mutual information
            mi = zeros(size(eeg_mi.trial{trial_idx}(cc,:)));
            sigAll = [];
            sigEff = [];
            % loop through all the peak times
            for pp = 1 : numel(ptimes_s)
                
                %current peak time
                tt = ptimes_s(pp);
                
                % the center is the sampling point that is closest to the
                % peak
                xcent = nearest(eeg_hg2.time{1}, tt);
                
                % concatenate all the data around this peak with a time
                % window that acommodates all lags and an additional
                % padding of half the window that is used to compute MI
                sigtmp = cat(1, ...
                    eeg_hg2.trial{1}(cc,...%high gamma
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2),...
                    eeg_lg2.trial{1}(cc,...%low gamma
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2),...
                    eeg_b2.trial{1}(cc,...%beta
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2),...
                    eeg_a2.trial{1}(cc,...%alpha
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2),...
                    eeg_t2.trial{1}(cc,...%theta
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2),...
                    eeg_d2.trial{1}(cc,...%delta
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2),...
                    eeg_2.trial{1}(cc,...%and raw
                    xcent-window_width/2-window_width_small/2:...
                    xcent+window_width/2+window_width_small/2));
                %concatenate all the peaks along the third dimension
                sigAll = cat(3, sigAll, sigtmp);
                
                %concatenate all the gamma signals on effect electrodes
                %along the third dimension only for the time window that is
                %used to compute MI
                sigEff = cat(3, sigEff,...
                    eeg_hgeff2.trial{1}(:,...
                    xcent-window_width_small/2:xcent+window_width_small/2));
                
            end
            
            % cut out the window for "zero lag" out of all data
            Ownshift = sigAll(:, window_width/2+1:...
                window_width/2 + window_width_small+1,:);
          
            warning off
            % now loop through all shifts from -window_width/2 to
            % +window_width/2
            for ssh = 0:1:window_width
             
                % the shifted signal of length "window_width_small"
                Sigshift = sigAll(:,1+ssh:ssh+window_width_small+1,:);
                
                %
                try
                    %compute the MI between shifted Signal and Effect
                    %signal conditioned on the "zero lag" signal
                    mi(1+ssh) = gccmi_ccc(Sigshift(:,:)', sigEff(:,:)', Ownshift(:,:)');
                catch
                    % error if signal is also "zero lag"
                    mi(1+ssh) = nan;
       
                end
            end
            
            % store this signal
            eeg_mi.trial{trial_idx}(cc,:) = mi;
        end
 
    end
    
    save(fullfile(deriv_pth3, ['LaggedMIaround_predictionPeaksFinal1s_slidingWindow.mat']), 'eeg_mi');
    
    
end
