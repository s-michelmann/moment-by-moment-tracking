%% load in the cut raw data.
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
%% Prepare the audio once for everyone!

interpolate_gamma_spikes = true;
order_fixed = false;
morder = nan;
mo_max = 50; %usually 50 ends on the same with higher max

audio_filename = [data_pth_ filesep 'stimuli' filesep 'Story_Original_MRI.wav'];

% read in the data
[audio_data, audio_fq] = audioread(audio_filename);
ft_sound = struct();

% create a fieldtrip structure for easy data handling
ft_sound.fsample = audio_fq;
ft_sound.time{1} =  [0:1/audio_fq:length(audio_data)/audio_fq-1/audio_fq];
ft_sound.trial{1} = mean(audio_data,2)';
ft_sound.label = {'audio'};
ft_sound.sampleinfo = [1, length(audio_data)];

%bpfilter between 200 and 5000Hz
cfg = [];
cfg.bpfilter = 'yes';
cfg.bpfreq = [200 5000];
audio_bp = ft_preprocessing(cfg, ft_sound);

%take the absolute of the hilbert transform and demean (detrend doesn't
%make sense for the clean digital audio!)
audio_ds = audio_bp;
audio_ds.trial{1} = abs(hilbert(ft_sound.trial{1}));
audio_ds.trial{1} = audio_ds.trial{1} - mean(audio_ds.trial{1});

%downsample to 100Hz to match EEG
cfg = [];
cfg.resamplefs = 100;
audio_ds = ft_resampledata(cfg, audio_ds);

% already define signal 3 here!
s3 = audio_ds.trial{1};

for ss = 1 : size(sjs, 1)
    
    %% prepare the EEG and the audio
    %subject code
    code_ = sjs(ss,:);
    
    % create correct path-names
    subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
    deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
        'sub-' code_ ];
    load(fullfile(deriv_pth, 'eeg_manualica_notch.mat'));
    deriv_pth = [data_pth_ 'derivatives' filesep 'audio_to_brain', filesep ...
        'sub-' code_ ];
    deriv_pth2 = [data_pth_ 'derivatives' filesep 'granger_final', filesep ...
        'sub-' code_ ];
    if ~exist(deriv_pth,'dir'); mkdir(deriv_pth); end
    if ~exist(deriv_pth2,'dir'); mkdir(deriv_pth2); end
    
    % select the notch filtered preprocessed post-ica EEG as eeg
    eeg = eeg_manualica;
    clear eeg_manualica;
    
    % bpfilter eeg in the gamma band (70-200Hz) note that 180 bs comes into
    % play here
    cfg = [];
    cfg.bpfilter = 'yes';
    cfg.bpfreq = [70 200];
    eeg_bp = ft_preprocessing(cfg, eeg);
    
    % loop through all the channels
    for cc = 1 : numel(eeg_bp.label)
        % and replace the data with the amplitude of the hilbert
        eeg_bp.trial{1}(cc, :) = abs(hilbert( eeg_bp.trial{1}(cc, :)));
        eeg_bp.trial{2}(cc, :) = abs(hilbert( eeg_bp.trial{2}(cc, :)));
    end
    
    % redefine the trial to ensure no trailing nans and near equal length
    % to audio
    cfg = [];
    cfg.toilim = [ft_sound.time{1}(1) ft_sound.time{1}(end)];
    eeg_bp = ft_redefinetrial(cfg, eeg_bp);
    
    % now downsample to 100Hz
    cfg = [];
    cfg.resamplefs = 100;
    eeg_bp = ft_resampledata(cfg, eeg_bp);
    
    
    % finally detrend the data to remove drift (but mostly mean!)
    cfg = [];
    cfg.detrend = 'yes';
    eeg_bp = ft_preprocessing(cfg, eeg_bp);
    
    
    %% Parameters FOR GC analysis
    
    regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
    icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)
    momax     = mo_max;     % maximum model order for model order is 50, which corresponds to 500ms
    
    
    %define datastructures to collect the results
    model_eeg = [];
    model_eeg.fsample = 100;
    
    model_eeg.trial{1} = zeros(numel(eeg_bp.label), length(s3));
    model_eeg.time{1} = eeg_bp.time{1};
    model_eeg.trial{2} = zeros(numel(eeg_bp.label), length(s3));
    model_eeg.time{2} = eeg_bp.time{1};
    model_eeg.trial{3} = zeros(numel(eeg_bp.label), length(s3));
    model_eeg.time{3} = eeg_bp.time{1};
    
    model_eeg.trial{4} = zeros(numel(eeg_bp.label), length(s3));
    model_eeg.time{4} = eeg_bp.time{1};
    model_eeg.trial{5} = zeros(numel(eeg_bp.label), length(s3));
    model_eeg.time{5} = eeg_bp.time{1};
    model_eeg.trial{6} = zeros(numel(eeg_bp.label), length(s3));
    model_eeg.time{6} = eeg_bp.time{1};
    model_eeg.label = eeg_bp.label;
    
    model_eeg.A = cell(numel(eeg_bp.label), 1);
    model_eeg.F = cell(numel(eeg_bp.label), 1);
    model_eeg.effect = cell(numel(eeg_bp.label), 1);
    
    % loop through all the channels
    for cc = 1 : numel(eeg_bp.label)
        
                
        fprintf('\nsj = %d channel = %d\n',ss, cc);
        % define signal 1 as the first trial
        s1 =  (eeg_bp.trial{1}(cc, 1:length(s3)));% the length should be the same anyway...
        
        
        
       
        
        % define signal 2 as the second trial
        s2 = (eeg_bp.trial{2}(cc, 1:length(s3)));% the length should be the same anyway...
        
        if interpolate_gamma_spikes % this fixes problems with residual electrical spikes
            s1 = despike_dilate(s1', 5, 5, false, 0)';
            s2 = despike_dilate(s2', 5, 5, false, 0)';
        end
        
		%save the signal 1 in trial 1
        model_eeg.trial{1}(cc,:) = s1;
		
        %save the signal 2 in trial 2
        model_eeg.trial{2}(cc,:) = s2;
        
        %save the audio in trial 3
        model_eeg.trial{3}(cc,:) = s3;
        
        
        %concatenate the data along the first dimension in X
        X = [s1;s2; s3];
        
        %% Model order estimation
        
        % Calculate information criteria up to max model order
        if ~order_fixed
            ptic('\n*** tsdata_to_infocrit\n');
            [AIC,BIC] = tsdata_to_infocrit(X,momax,icregmode, false);
            ptoc('*** tsdata_to_infocrit took ');
            
            [~,bmo_AIC] = min(AIC); % use Akaice Information Criterion for model selection
            fprintf('\nbest model order (AIC) = %d\n',bmo_AIC);
            
            % Select model order
            morder = bmo_AIC;
        end
        ptic('\n*** GCCA_tsdata_to_pwcgc... ');
               
        % finally estimate Granger causality in the traditional way i.e in
        % the time domain
        [F,A,SIG, E] = GCCA_tsdata_to_pwcgc(X,morder,regmode); % use same model order for reduced as for full regressions
        ptoc;
        
        % save the residuals of signal 1 in trial 4
        model_eeg.trial{4}(cc,morder+1:end) = E(1,:);
        
        % save the residuals of signal 2 in trial 5
        model_eeg.trial{5}(cc,morder+1:end) = E(2,:);
        
        % save the residuals of the audio in trial 6
        model_eeg.trial{6}(cc,morder+1:end) = E(3,:);
        
        model_eeg.A{cc} = A;
        model_eeg.F{cc} = F;
        
        %the first index |i| of |F| is the target (causee) variable, ...
        % the second |j| the source (causal) variable
        % F(1,2) is s2 'causing' s1, f(2,1) is s1 'causing' s2
        model_eeg.effect{cc} = F(1,2) - F(2,1);
        
        % Check for failed (full) regression
        assert(~isbad(A),'VAR estimation failed');
        
        % Check for failed GC calculation
        assert(~isbad(F,false),'GC calculation failed');
        
        % Check VAR parameters (but don't bail out on error - GCCA mode is quite forgiving!)
        rho = var_specrad(A);
        fprintf('\nspectral radius = %f\n',rho);
        if rho >= 1,       fprintf(2,'WARNING: unstable VAR (unit root)\n'); end
        if ~isposdef(SIG), fprintf(2,'WARNING: residuals covariance matrix not positive-definite\n'); end
              
    end
    % note: other versions e.g. without inclusion of the audio are 
    % essentially the same result.
    % spike interpolation is only absolutely necessary for a few 
    % (almost dead electrodes) in a single subject. Those electrodes 
    % could also be excluded but interpolation improves overall result, so 
    % that's preferred. Other frequencies don't show results, using raw 
    % data seems to be a weaker version of the gamma result.
    save(fullfile(deriv_pth2, ['CrossRunInclAudio_gamma70_200Hz_ManualICAnotchMO50max_interpGammaSpike.mat']), 'model_eeg');
    
end

