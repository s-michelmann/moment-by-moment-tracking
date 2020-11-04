clear;
clc

% 205 subjects were recorded
% for the brain correlation
% 13 subjects were rejected because they had a cosine-distance of >.85 to
% the average of all other subjects in one of the runs
load('dataBehavioralSegmentation.mat')
%% collect the data in a vector of ones
close all
dts1 = zeros(size(datcell,1), 7.5*60000);
dts2 = zeros(size(datcell,1), 7.5*60000);
% figure;
zero_vec = zeros(7.5*60000, 1);

% data are integrated within a 1 second window...
sm_ = 1000;

for s = 1:size(datcell,1) 
    s
    % write ones in a vector of zeros
    zero_vec = zeros(7.5*60000, 1);
    zero_vec(datcell{s, 1}) = 1;
    %smooth and threshold
    zero_vec = smooth(zero_vec, sm_);
    zero_vec(zero_vec>0) = 1;
    
    % subject contributes ones in the second around the button press
    dts1(s,:) = (zero_vec);
    
    % the same for the second run
     zero_vec = zeros(7.5*60000, 1);
     zero_vec(datcell{s, 2}) = 1;
     zero_vec = smooth(zero_vec, sm_);
     zero_vec(zero_vec>0) = 1;
     dts2(s,:) = (zero_vec);
     
end
%% convert to cosine similarity instead!!! (for consistency)
% do this later after the randomization has finished!;
%% QQ plot
close all

line([prctile(mean(dts1)', 95) prctile(mean(dts1)', 95)],[-.04 0.44], 'color', [40 156 174]./255,'LineStyle','-');
hold on
line([prctile(mean(dts1)', 99) prctile(mean(dts1)', 99)],[-.04 0.44], 'color', [40 156 174]./255,'LineStyle','-');

line([prctile(mean(dts1)', 97.5) prctile(mean(dts1)', 97.5)],[-.04 0.44], 'color', [40 156 174]./255,'LineStyle','-');
line([prctile(mean(dts1)', 99.9) prctile(mean(dts1)', 99.9)],[-.04 0.44], 'color', [40 156 174]./255,'LineStyle','-');
set(gca,'fontsize', 12);
% qqplot(mean(dts1), mean(dts2)); xlim([-.04 0.44])
h = qqplot(mean(dts1), mean(dts2),[1:1/1000:100]); xlim([-.04 0.44]);ylim([-.04 0.44])
set(h(1), 'Marker', 'o', 'markersize',4,'markeredgecolor',[116 78 144]./255 );
set(h, 'color',[0 0 0]);
xlabel('run 1 agreement');
ylabel('run 2 agreement');
set(gcf, 'color', 'w');
%

%% now check the distributions (come back to kurtosis later)

signal1 = mean(dts1)';
signal2 = mean(dts2)';
n = length(signal1);
dtsall = [dts1; dts2];
sig_all = mean(dtsall);

maxmin_range = max(sig_all)-min(sig_all);

%get binsize via FD rule
fd_bins      = ceil(maxmin_range/(2.0*iqr(sig_all)*n^(-1/3))); % Freedman-Diaconis 

close all;
clc
figure
ha = tight_subplot(1, 2, 0);
axes(ha(1));
[counts,bins] = hist(mean(dts1)', fd_bins); %# get counts and bin locations
h = barh(bins,counts, 'FaceColor', [30 83 145]./255)
set(get(h,'Parent'),'xdir','r');
ylim([0 0.35])
ylabel('agreement');

axes(ha(2));
ylim([0 0.35]);

[counts,bins] = hist(mean(dts2)', fd_bins); %# get counts and bin locations
h = barh(bins,counts, 'FaceColor', [198 0 27]./255);
set(gca,'ytick',[]);
set(gcf, 'color', 'w')

%% do the permutation

signal1 = mean(dts1)';
signal2 = mean(dts2)';
n = length(signal1);
dtsall = [dts1; dts2];
sig_all = mean(dtsall);
maxmin_range = max(sig_all)-min(sig_all);
%get binsize via FD rule
fd_bins      = ceil(maxmin_range/(2.0*iqr(sig_all)*n^(-1/3))); % Freedman-Diaconis 

% bin data (via Matlab function hist)
nbins = fd_bins;
[hdat1,x1] = hist(signal1,nbins);
[hdat2,x2] = hist(signal2,nbins);

% convert histograms to probability values
hdat1 = hdat1./sum(hdat1);
hdat2 = hdat2./sum(hdat2);

run1k = kurtosis(hdat1);
run2k = kurtosis(hdat2);

% also check entropy? 
run1E = -sum(hdat1.*log2(hdat1+eps));
run2E = -sum(hdat2.*log2(hdat2+eps));

n_rand = 1000;

avg_mat = [1/size(dts2,1).*ones(1, size(dts2,1)),...
    zeros(1, size(dts2,1));...
    zeros(1, size(dts2,1)) ...
    1/size(dts2,1).*ones(1, size(dts2,1))];
run1Er    = zeros(n_rand, 1);
run2Er    = zeros(n_rand, 1);
rundiffEr = zeros(n_rand, 1);
run1kr   = zeros(n_rand, 1);
run2kr    = zeros(n_rand, 1);
rundiffkr = zeros(n_rand, 1);

for rr = 1 : n_rand
    rr
    prm_ = randperm(size(dtsall, 1))'; 
    pr_mat = eye(size(dtsall,1));
    pr_mat = pr_mat(prm_,:);
    sig12 = avg_mat*pr_mat*dtsall;   
    [hdatr1,x1] = hist(sig12(1,:),nbins);
    [hdatr2,x2] = hist(sig12(2,:),nbins);    
    % convert histograms to probability values
    hdatr1 = hdatr1./sum(hdatr1);
    hdatr2 = hdatr2./sum(hdatr2);    
    run1Er(rr) = -sum(hdatr1.*log2(hdatr1+eps));
    run2Er(rr) = -sum(hdatr2.*log2(hdatr2+eps));
    rundiffEr(rr) = run2Er(rr)-run1Er(rr);
    run1kr(rr) = kurtosis(hdatr1);
    run2kr(rr) = kurtosis(hdatr2);
    rundiffkr(rr) = run2kr(rr)-run1kr(rr);
end

%%
p_val = numel(find(rundiffkr > (run2k - run1k)))/n_rand;



figure;
Violin(rundiffkr, 1, 'ViolinColor', [207 218 230]./255, 'Width', 0.3);
hold on;
plot(0.75, (run2k - run1k), 'o', 'MarkerEdgeColor',  [27 0 198]./255,...
    'MarkerFaceColor', [116 78 144]./255);
xlim([-1 3])
ylim([-15 15])
xlabel('random permutation of run 2 and 1')
ylabel('difference in kurtosis')
set(gcf, 'color', 'white');
set(gca,'fontsize', 12);

%%
%lool compute the distances for each subject
dist11 = zeros(length(datcell),1);
% dist12 = zeros(length(datcell),1);
% dist21 = zeros(length(datcell),1);
dist22 = zeros(length(datcell),1);
for ss = 1 : length(datcell)
    avg_vec = ones(1, length(datcell))./(length(datcell)-1);
    ss
    avg_vec(ss) = 0;
    dist_1 = avg_vec* dts1;  
    dist_2 = avg_vec* dts2;   
    
    vec_1 = dts1(ss, :); 
    vec_2 = dts2(ss, :); 
    
     dist11(ss) = pdist2(dist_1, vec_1, 'cosine');
%     dist12(ss) = pdist2(dist_1, vec_2, 'cosine');
%     dist21(ss) = pdist2(dist_2, vec_1, 'cosine');
    dist22(ss) = pdist2(dist_2, vec_2, 'cosine');
end

mat(1,1) = mean(dist11);
% mat(2,1) = mean(dist21);
% mat(1,2) = mean(dist12);
mat(2,2) = mean(dist22);

%%
plot([1 20] , [dist11 dist22], '--o', 'color', [.6 .6 .6]);
xlim([0 21])
hold on;
plot([1 20] , [mean(dist11) mean(dist22)], '-o', 'color', [1 0 0], 'LineWidth', 2);
%% plot similarity instead (1- dist)
figure
plot([0 1], [0 1]); hold on
plot(1-dist11, 1-dist22, 'ko');
hold on;
plot(mean(1-dist11), mean(1-dist22), 'r*', 'MarkerSize', 15);
set(gcf, 'color', 'w');
xlabel('cosine similarity to others run 1')
ylabel('cosine similarity to others run 2')
% %
xlim([0 0.65])
ylim([0 0.65])
%%
% imagesc(mat)
% 
agreement_change = mean(dist11-dist22);
% self_agreement_change = mean(dist12-dist22); or similar

numel(find(dist22 < dist11))/205

% probability to get that many ppl with .5 chance of improvement
1 - binocdf(numel(find(dist22 < dist11)),205,0.5)
%% random permutation of ISC

n_rand = 1000;
N = length(datcell);
distdifr = zeros(n_rand,1);

przlearners = zeros(n_rand,1);


for rr = 1 : n_rand
    tic;
    rr
  
    prm_ = randperm(size(dtsall, 1))';
    dist11r = zeros(N,1);
    dist22r = zeros(N,1);
    
    
    pr_mat = eye(size(dtsall,1));
    pr_mat = pr_mat(prm_,:);
    
    dt_pr = pr_mat*dtsall;
%     sel1 = pr_mat(1:size(dts1, 1),:);
%     sel2 = pr_mat(size(dts1, 1)+1:end,:);
    
    avg_vec1 = [ones(1, N)./(N-1), zeros(1,N)];
    avg_vec2 = [zeros(1,N), ones(1, N)./(N-1)];
    for ss = 1 : N
        
        
        avg_vec1(ss) = 0;
        avg_vec2(ss+N) = 0;
        dist_1 = avg_vec1* dt_pr;% mean(dts1([1:ss-1, ss+1:end], :));
        dist_2 = avg_vec2* dt_pr;% m
        
        vec_1 = dt_pr(ss,:) ;%sel1(ss,:)* dtsall; 
        vec_2 = dt_pr(ss+N,:) ;%sel2(ss,:)* dtsall; %% question mark!!
        
        dist11r(ss) = pdist2(dist_1, vec_1, 'cosine');
        dist22r(ss) = pdist2(dist_2, vec_2, 'cosine');
        avg_vec1(ss) = 1/(N-1);
        avg_vec2(ss+N) = 1/(N-1);
    end
    
    distdifr(rr) = mean(dist11r-dist22r);
    przlearners(rr) = numel(find(dist22r < dist11r))/N;
    toc
end

save distdifr distdifr
save przlearners przlearners

%%

figure;
mean(dist11-dist22)

numel(find(dist22 < dist11))/N


p_val = numel(find(distdifr > mean(dist11-dist22)))/n_rand;



figure;
Violin(distdifr, 1, 'ViolinColor', [40 156 174]./255, 'Width', 0.3);
hold on;
plot(0.75, ( mean(dist11-dist22)), 'o', 'MarkerEdgeColor',  [27 0 198]./255,...
    'MarkerFaceColor', [40 156 174]./255);


xlim([-1 3])
ylim([-0.06 0.06])
xlabel('random permutation of run 2 and 1')
ylabel('increase in cosine similarity to others ')
set(gcf, 'color', 'white');
set(gca,'fontsize', 12);


nhist(distdifr, 'color', [151 243 255 ]./255);
set(gcf, 'color', 'w');
line([mean(dist11-dist22) mean(dist11-dist22)],[0 120], 'color', 'red','LineStyle','--', 'LineWidth',2);
xlabel('Change in cosine distance to others');
set(gca,'fontsize', 12);
% xlim([-20 20])
%% 


mean(dist11-dist22)

prz_real = numel(find(dist22 < dist11))/N


p_val2 = numel(find(przlearners > prz_real))/n_rand;
nhist(przlearners, 'color', [151 243 255 ]./255);
set(gcf, 'color', 'w');
hold on
line([prz_real prz_real],[0 120], 'color', 'red','LineStyle','--', 'LineWidth',2);
xlabel('ratio learners');
set(gca,'fontsize', 12);
xlim([0.35 0.65])