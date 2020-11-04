function Main_Server_MI2D_Event_Boundary_connectivity(nr, trial_idx)

sjs = [''; ''; '';''; ''; ''; ''; ''; '']; % codes removed for sharing
code_ = sjs(nr,:);



load('ptimes_s.mat');
load('all_models.mat');
addpath(genpath('gcmi-master'));
addpath('fieldtrip');
ft_defaults;
% cd(code_);
load(fullfile(code_, 'eeg_aall.mat'));
load(fullfile(code_, 'eeg_all.mat'));
load(fullfile(code_, 'eeg_ball.mat'));
load(fullfile(code_, 'eeg_dall.mat'));
load(fullfile(code_, 'eeg_hgall.mat'));
load(fullfile(code_, 'eeg_hgeff.mat'));
load(fullfile(code_, 'eeg_lgall.mat'));
load(fullfile(code_, 'eeg_mi1.mat'));
load(fullfile(code_, 'eeg_mi2.mat'));
load(fullfile(code_, 'eeg_tall.mat'));

fid = fopen([code_ filesep 'HC_channels.tsv'], 'r');
txt = textscan(fid, '%s', 'delimiter', '\t', 'Headerlines', 0);
fclose all   ;
%% now compute the mutual information at peak-locked data
% define the window of possible lags (-2.5s to +2.5s)
window_width_lag = 5*eeg_hgeff.fsample;
window_width_shift = 6*eeg_hgeff.fsample;
% define the window with to compute mutual information across
window_width_small = eeg_hgeff.fsample;

sampling_rate = all_models{1}.fsample;

try
eeg_mi2  = rmfield(eeg_mi2, 'sampleinfo');
catch
end
try
eeg_mi1  = rmfield(eeg_mi1, 'sampleinfo');
catch
end
for tt = 1 : window_width_shift + 1
    eeg_mi1.trial{tt} = zeros(numel(txt{1}), ...
        window_width_lag+1);
    eeg_mi1.time{tt} = [-window_width_lag/2:100/eeg_all.fsample:window_width_lag/2];
    eeg_mi2.trial{tt} = zeros(numel(txt{1}), ...
        window_width_lag+1);
    eeg_mi2.time{tt} = [-window_width_lag/2:100/eeg_all.fsample:window_width_lag/2];

end
eeg_mi1.label = txt{1};
eeg_mi2.label = txt{1};
rng(rand_nr);
%
% % resize the boundary time course to the EEG SR
% tc2 = imresize(rand_tc, [1 45000]);
% % select only the times from 0(interp = 0.002?) to end
% tc2 = tc2(:, 1:length(all_models{nr}.time{1}));
% padding to exclude peaks around the bounds!
pad_sec = 3;


time_pick = true(size(all_models{nr}.time{1}));
time_pick(1:(sampling_rate*pad_sec)) = 0;
time_pick(end-(sampling_rate*pad_sec-1):end) = 0;

all_time_s =  all_models{nr}.time{1}(time_pick);


eeg_hg2.trial =  eeg_hgall.trial(trial_idx); % = ft_selectdata(cfg, eeg_hgall);
eeg_lg2.trial =  eeg_lgall.trial(trial_idx); % = = ft_selectdata(cfg, eeg_lgall);
eeg_b2.trial =  eeg_ball.trial(trial_idx); % = = ft_selectdata(cfg, eeg_ball);
eeg_a2.trial =  eeg_aall.trial(trial_idx); % = = ft_selectdata(cfg, eeg_aall);
eeg_t2.trial =  eeg_tall.trial(trial_idx); % = = ft_selectdata(cfg, eeg_tall);
eeg_d2.trial =  eeg_dall.trial(trial_idx); % = = ft_selectdata(cfg, eeg_dall);
eeg_2.trial =  eeg_all.trial(trial_idx); % =  = ft_selectdata(cfg, eeg_all);
eeg_hgeff2.trial =  eeg_hgeff.trial(trial_idx); % = = ft_selectdata(cfg, eeg_hgeff);

eeg_hg2.time =  eeg_hgall.time(trial_idx); % = ft_selectdata(cfg, eeg_hgall);
eeg_lg2.time =  eeg_lgall.time(trial_idx); % = = ft_selectdata(cfg, eeg_lgall);
eeg_b2.time =  eeg_ball.time(trial_idx); % = = ft_selectdata(cfg, eeg_ball);
eeg_a2.time =  eeg_aall.time(trial_idx); % = = ft_selectdata(cfg, eeg_aall);
eeg_t2.time =  eeg_tall.time(trial_idx); % = = ft_selectdata(cfg, eeg_tall);
eeg_d2.time =  eeg_dall.time(trial_idx); % = = ft_selectdata(cfg, eeg_dall);
eeg_2.time =  eeg_all.time(trial_idx); % =  = ft_selectdata(cfg, eeg_all);
eeg_hgeff2.time =  eeg_hgeff.time(trial_idx); % = = ft_selectdata(cfg, eeg_hgeff);
   
delete(gcp('nocreate'))

parpool('local',16)
cc_new = 1;
% loop through channels
for cc = 1 : numel(eeg_hgall.label)

    tic;
    % this is for restriction to the HC
%     if ismember(eeg_hgall.label{cc}, txt{1})
%         fprintf('\nsj = %d trial = %d channel = %d',ss, trial_idx, cc);
%         
%     else
%         continue;
%     end


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
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2),...
            eeg_lg2.trial{1}(cc,...%low gamma
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2),...
            eeg_b2.trial{1}(cc,...%beta
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2),...
            eeg_a2.trial{1}(cc,...%alpha
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2),...
            eeg_t2.trial{1}(cc,...%theta
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2),...
            eeg_d2.trial{1}(cc,...%delta
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2),...
            eeg_2.trial{1}(cc,...%and raw
            xcent-window_width_lag/2-window_width_small/2-window_width_shift/2:...
            xcent+window_width_lag/2+window_width_small/2+window_width_shift/2));
        %concatenate all the peaks along the third dimension
        sigAll = cat(3, sigAll, sigtmp);

        %concatenate all the gamma signals on effect electrodes
        %along the third dimension only for the time window that is
        %used to compute MI
        sigEff = cat(3, sigEff,...
            eeg_hgeff2.trial{1}(:,...
            xcent-window_width_small/2-window_width_shift/2:...
            xcent+window_width_small/2+window_width_shift/2));

    end

    ln_ = length(eeg_mi1.trial{trial_idx}(1,:));
    mis_ = zeros(window_width_shift+1, length(eeg_mi1.trial{trial_idx}(1,:)));
    % loop that goes through all points around event boundary
    parfor llg = 0 : window_width_shift

        % cut out the window for "zero lag" out of all data
        % own shift should always be at the same lag as Effect
        Ownshift = sigAll(:, llg + window_width_lag/2+1:...
            llg + window_width_lag/2+ window_width_small+1,:);
        Efflag = sigEff(:,llg+1:llg+window_width_small+1,:);
        warning off

        % now loop through all lags from -window_width/2 to
        % +window_width/2
        %initialize the mutual information
        mi = zeros(1, ln_);

        for ssh = 0:1:window_width_lag
            %%  ssh
            % the shifted signal of length "window_width_small"
            Sigshift = sigAll(:,1+ssh+llg:ssh+llg+window_width_small+1,:);

            %
            try
                %compute the MI between shifted Signal and Effect
                %signal conditioned on the "zero lag" signal
                mi(1+ssh) = gccmi_ccc(Sigshift(:,:)', Efflag(:,:)', Ownshift(:,:)');
            catch
                % error if signal is also "zero lag"
                mi(1+ssh) = nan;

            end

        end
       mis_(llg+1,:) = mi;

    end
%     toc
    if trial_idx== 1
        for llg = 0 : window_width_shift
            eeg_mi1.trial{llg+1}(cc_new,:) = mis_(llg+1,:) ;
        end
        cc_new = cc_new+1;
    else
         for llg = 0 : window_width_shift
            eeg_mi2.trial{llg+1}(cc_new,:) = mis_(llg+1,:) ;
         end
         cc_new = cc_new+1;
    end
toc
end


if trial_idx == 1
    save(fullfile(code_, ['2DMIaround_Boundaries'  num2str(rand_nr) '.mat']), 'eeg_mi1', '-v7.3');

else
    save(fullfile(code_, ['2DMIaround_Boundaries'  num2str(rand_nr) '.mat']), 'eeg_mi2', '-v7.3');

end
end
