%% note: codes removed for sharing, lots of plotting and manual checking 
% + sanity checks removed for structure and readability 
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
code_ = sjs(1,:);
subj_pth_ = [data_pth_ 'sub-' code_ filesep 'ieeg'];
% define trials
deriv_pth = [data_pth_ 'derivatives' filesep 'preprocessing', filesep ...
    'sub-' code_ ];

load([deriv_pth filesep 'eeg_cut.mat'])

%% read in electrodes with known coordinates
fid = fopen([subj_pth_ filesep 'sub-' code_ '_electrodes.tsv']);
txt = textscan(fid, '%s%f%f%f%s', 'delimiter', '\t', 'Headerlines', 1);
fclose all;

%% despike to deal with artifacts

% 0 automatic artifact/spike detection based on criteria: 
    %>5iq range >3 channels spike together in a +- 300ms range
    
% 1st identify artifactual areas and exclude them for ICA (correct
% manually)
% 2nd do ICA reref with a conservative (manual) exclusion
% 3rd interpolate spikes only if they fall within the manual artifact
% areas (+-2sec) with a liberal threshold of 3.5 IQ and padding around it

% sampling info in raw file corresponding to data
si_vector = [eeg.sampleinfo(1,1):eeg.sampleinfo(1,2), ...
    eeg.sampleinfo(2,1):eeg.sampleinfo(2,2)];

% after some semi-automatic pre-checking, do the actual artifact definition 
% manually in the fieldtrip data browser!
cfg.blocksize = 25;
cfg.viewmode = 'vertical';
cfg.detrend = 'yes';
cfg.continuous = 'yes';
cfg_art_range = ft_databrowser(cfg, tmpeeg); 
%%
%cfg_art_range has the information about artifactual samples
save([deriv_pth filesep 'cfg_art_range.mat'], 'cfg_art_range'); 
%%
%% do the filter for later ICA (we don't want line noise to pull components)
tmpeeg = eeg;
cfg = [];
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [55 65; 115 125; 175 185]; % this is because noise is too strong here... (refilter later bc none of those will be done to the data)
tmpeeg             = ft_preprocessing(cfg, tmpeeg);
%%
data_vec = horzcat(tmpeeg.trial{:})';
art_vec = zeros(size(si_vector));
for aa = 1: size(cfg_art_range.artfctdef.visual.artifact,1)
    
    art_vec(cfg_art_range.artfctdef.visual.artifact(aa,1)-0.2*eeg.fsample <= ...
        si_vector & cfg_art_range.artfctdef.visual.artifact(aa,2)+0.2*eeg.fsample >= ...
        si_vector) = 1;
end

% overwrite with nan for ICA computation
data_vec(find(art_vec), :) = nan;
save([deriv_pth filesep 'art_vec.mat'], 'art_vec');

%% do the ica 

data_comp       = ft_componentanalysis([], tmpeeg); %using eeglab runica

unmixing        = data_comp.unmixing;
used_labels     = eeg.label;
mixing          = inv(unmixing);

save([deriv_pth filesep 'unmixing.mat'], 'unmixing');
save([deriv_pth filesep 'used_labels.mat'], 'used_labels');
% mixing * components = data
% unmixing * data = components 

cfg          = [];
cfg.viewmode = 'vertical';
cfg.blocksize = 25;
cfg.continuous = 'yes';
ft_databrowser(cfg,data_comp);

% check the distribution of hidden sources;
figure;
imagesc(abs(mixing))
caxis([-max(abs(mixing(:))) max(abs(mixing(:)))]);
colormap(jet(256)); 
colorbar;

set(gca, 'XTick',1:numel(used_labels));
set(gca, 'YTick',1:numel(used_labels));
set(gca, 'YTickLabel', used_labels);
set(gcf, 'Color', 'w');
%% 
close all
%% interpolate the cut data

load([deriv_pth filesep 'art_vec.mat']);
load([deriv_pth filesep 'eeg_cut.mat']);
% w +-150ms padding

cfg = [];
cfg.preproc.detrend = 'yes'; 
eeg_interp = ft_preprocessing(cfg, eeg);

exclusion_range = round(0.15*eeg.fsample);

data_vec = horzcat(eeg.trial{:})';

%despike_dilate(x, spk_thresh, exclusion_range, do_log, spk_side, permission_range)
[data_despike, spikes] = despike_dilate(data_vec,3.5, exclusion_range, 0, 0, art_vec');
% data_despike(spikes) = nan;

eeg_interp.trial{1} = data_despike(1:size(eeg.trial{1},2), :)';
eeg_interp.trial{2} = data_despike(size(eeg.trial{1},2)+1:end, :)';
%%

cfg = [];
cfg.blocksize = 25;
cfg.viewmode = 'vertical';
cfg.detrend = 'yes';
cfg.continuous = 'yes';
ft_databrowser(cfg, eeg_interp);
%%
save([deriv_pth filesep 'eeg_interp_raw.mat'], 'eeg_interp');

%% now apply the unmixing matrix to he original eeg 
load([deriv_pth filesep 'unmixing.mat']);
mixing = inv(unmixing);
used_labels = eeg.label;
cfg             = [];
cfg.topolabel   = used_labels;
cfg.unmixing    = unmixing;
data_comp       = ft_componentanalysis(cfg, eeg_interp);
%%
close all

chi_val = sum(bsxfun(@rdivide, ((bsxfun(@minus,abs(mixing),  mean(abs(mixing),1))).^2), ...
    mean(abs(mixing),1)),1);

%sort the chi-square values and keep track of the sorting-indices
[chisort, s_indices] = sort(chi_val, 'descend');


%% plot again the mixing matrix(topography) with components sorted by broadness
figure
imagesc((mixing(:, s_indices)))
caxis([-max(abs(mixing(:))) max(abs(mixing(:)))]); %axis xy
colormap(jet(256)); 
colorbar;

set(gca, 'XTick',1:numel(used_labels))
set(gca, 'YTick',1:numel(used_labels))

set(gca, 'YTickLabel', used_labels)
tmp = [1:numel(used_labels)];
set(gca, 'XTickLabel', num2cell(tmp(s_indices)))
colormap(jet(256));

data_comp.trial{1}(:, any(isnan(data_comp.trial{1}))) = 0;
data_comp.trial{2}(:, any(isnan(data_comp.trial{2}))) = 0;
%% do a manual threshold
reject_components   = find(~ismember(1:size(data_comp.label), s_indices(1:end-2))); %reject 0-3 broad components
cfg.component       = reject_components;
eeg_manualica             = ft_rejectcomponent(cfg, data_comp);
%% save
save([deriv_pth filesep 'eeg_manualica.mat'], 'eeg_manualica');
%% look
cfg = [];
cfg.blocksize = 25;
cfg.viewmode = 'vertical';
cfg.detrend = 'yes';
cfg.continuous = 'yes';
ft_databrowser(cfg, eeg_manualica);

%% Save a denotched version of this dataset
cfg = [];
cfg.bsfilter    = 'yes';
cfg.bsfreq      = [58 62; 118 122; 178 182]; % this is because noise is too strong here...
eeg_manualica = ft_preprocessing(cfg, eeg_manualica);

cfg = [];
cfg.blocksize = 25;
cfg.viewmode = 'vertical';
cfg.detrend = 'yes';
cfg.continuous = 'yes';
ft_databrowser(cfg, eeg_manualica);

%% Save a denotched version of this dataset
save([deriv_pth filesep 'eeg_manualica_notch.mat'], 'eeg_manualica');
