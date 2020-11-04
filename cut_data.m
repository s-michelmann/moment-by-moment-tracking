%% NOTE the electrodes.tsv is not in the same order as the data.label. This should be no problem if elec is selected based on labels!

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
%%

acq = ''; %'_acq-a'; %for some patients it's more complicated. this script is the standard case
code_ = sjs(1,:);
subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
    'sub-' code_ ];
if ~exist(deriv_pth,'dir'); mkdir(deriv_pth); end
% define trials

cfg            = [];
cfg.dataset    =  [subj_pth_ filesep 'sub-' code_ ...
    '_task-pieman' acq '_ieeg.edf'];

% read and preprocess the data
cfg.continuous = 'yes';
cfg.channel    = 'all';

% reading in the data
error = [];
try
    data           = ft_preprocessing(cfg);
catch error
    
    disp(error);
end

%% read in electrodes with known coordinates
fid = fopen([subj_pth_ filesep 'sub-' code_ '_electrodes.tsv']);
txt = textscan(fid, '%s%f%f%f%s', 'delimiter', '\t', 'Headerlines', 1)
fclose all;
%% select those channels as eeg
% the elec doesn't need to have the same order as the data.label!
data.elec.elecpos = horzcat(txt{2:4});
data.elec.label = txt{1};

cfg = [];
cfg.channel = txt{1};
eeg = ft_selectdata(cfg, data);

%% check
close all
ft_plot_mesh(eeg.elec.elecpos); 

%% check 2
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.detrend = 'yes';
ft_databrowser(cfg, eeg);

%% 
% print bad channels
% 
%  fid = fopen([deriv_pth filesep 'bad_channels.tsv'], 'w');

% fprintf(fid, '' );
% fprintf(fid, '\n');
% fprintf(fid, '');
% fprintf(fid, '\n');
%% OR load bad channels
fclose all;

fid = fopen([deriv_pth filesep 'bad_channels.tsv'], 'r');
txt = textscan(fid, '%s', 'delimiter', '\t', 'Headerlines', 0);
fclose all;

%% 
cfg = [];
cfg.channel = eeg.label(~ismember(eeg.label, ...
    txt{1}));
eeg = ft_selectdata(cfg, eeg);

%%
 
cfg = [];
cfg.channel = data.label(cellfun(@(x) contains(x, 'DC'), data.label));
trig = ft_selectdata(cfg, data);

%% indices are the first samplepoint where the flank of the trigger-channel 
% rises 
load([deriv_pth filesep 'trigger_indices.mat']);
%% sanity check (trigger onset and offset should have the same difference
disp(diff(indices')./(eeg.fsample*60))
%% now cut +-5 seconds
trl = [indices(:,1) - round(5  * eeg.fsample), ...
    indices(:,2) + round(5  * eeg.fsample), ...
    repmat(-round(5  * eeg.fsample), 2,1)];
cfg         = [];
cfg.trl = trl;
eeg = ft_redefinetrial(cfg, eeg);
%%
cfg = [];
cfg.viewmode = 'vertical';
cfg.preproc.detrend = 'yes';
ft_databrowser(cfg, eeg);

save([deriv_pth filesep 'eeg_cut.mat'], 'eeg');

