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


%% plot the little zoom-in inset
figure;
plot([1:450000]./60000,mean(dts1), 'color', [30 83 145]./255, 'LineWidth',2);
hold on;
plot([1:450000]./60000, mean(dts2), 'color', [198 0 27]./255, 'LineWidth',2);
hold on;
xlim([3 4]);
xlabel('time (minutes)');
ylabel('agreement');
set(gcf, 'color', 'white');
legend({'run 1',  'run 2'});
legend boxoff;
set(gca,'FontSize',20)

%%

close all
[C, LAG] = xcorr(mean(dts2),  mean(dts1), 'coeff');

[mm, ind] = max(C)
LAG(ind)
figure;
plot(LAG, C, 'k','LineWidth',2);
line([LAG(ind) LAG(ind)],[.2 1], 'color', 'red','LineStyle','--','LineWidth',2);
set(gca,'fontsize', 12);
xlim([-7000 7000])
ylim([0.4 1]);
set(gcf, 'color', 'w');

ylabel('correlation coefficient');
xlabel('lag milliseconds');
real_lag = LAG(ind)
%%
rand_N = 1000;
rand_lags = zeros(rand_N, 1);
dts_both = cat(1, dts1, dts2);
for rr = 1 : rand_N
    rr
    dts_both = dts_both(randperm(size(dts_both, 1))',:);
    dts2r = dts_both(1:size(dts2,1), :);
    dts1r = dts_both(size(dts1,1)+1 : end, :);
    [C, LAG] = xcorr(mean(dts2r),  mean(dts1r));
    % figure;
    % plot(LAG, C)
    [mm, ind] = max(C);
    rand_lags(rr) = LAG(ind);
end
%%

p_val = numel(find(rand_lags < real_lag))/rand_N;

figure;
Violin(rand_lags, 1, 'ViolinColor', [40 156 174]./255, 'Width', 0.3);
hold on;
plot(1, real_lag, 'o', 'MarkerEdgeColor',  [27 0 198]./255,...
    'MarkerFaceColor', [27 0 198]./255);
xlim([-1 3])
ylim([-200 200])
xlabel('random permutation of run 2 and 1')
ylabel('lag (ms)')
set(gcf, 'color', 'white');
set(gca,'fontsize', 12);
%%

%% make a selection based on the cosine distance to the mean of all others

%lool compute the distances for each subject
dist11 = zeros(length(datcell),1);
dist22 = zeros(length(datcell),1);
for ss = 1 : length(datcell)
    % this vector are the weights for a weighted sum (which is the average)
    avg_vec = ones(1, length(datcell))./(length(datcell)-1);
    ss
    % set this subjects weight to 0
    avg_vec(ss) = 0;
    
    % and multiply vector times matrix
    dist_1 = avg_vec* dts1;% mean(dts1([1:ss-1, ss+1:end], :));    
    dist_2 = avg_vec* dts2;   
    
    % now extract this subjects response vector
    vec_1 = dts1(ss, :); 
    vec_2 = dts2(ss, :); 
    
    % and compute the cosine similarity
    dist11(ss) = pdist2(dist_1, vec_1, 'cosine');
    dist22(ss) = pdist2(dist_2, vec_2, 'cosine');
end

%%
figure;
plot([1 20] , [dist11 dist22], '--o', 'color', [.6 .6 .6]);
xlim([0 21])
hold on;
plot([1 20] , [mean(dist11) mean(dist22)], '-o', 'color', [1 0 0], 'LineWidth', 2);
%%
cut_off = .85; %mean(distboth) + 2*std(distboth);
select_vec = (dist11 < cut_off) | (dist22 < cut_off);
reject_vec = (~select_vec);
%% 
sm_ = 2000;
figure;
smd1 = smoothdata(mean(dts1(select_vec,:)),'gaussian',sm_);
% plot(mean(dts1), 'k'); hold on;
plot(smd1, 'k--');
%
extract_events1 = bwlabeln(smd1 > prctile(smd1, 95));

hold on ;
smd1(extract_events1==0) = nan;
plot(smd1, 'r');
e1times = zeros(max(extract_events1),1);
for ee = 1: max(extract_events1)
    zvec = zeros(size(smd1));
    zvec(extract_events1 == ee) = smd1(extract_events1 == ee);
    [~, e1times(ee)] = max(zvec);
end
max(extract_events1)
%% 
figure;
smd2 = smoothdata(mean(dts2(select_vec,:)),'gaussian',sm_);
% plot(mean(dts2), 'k'); hold on;
plot(smd2, 'k--');
%
extract_events2 = bwlabeln(smd2 > prctile(smd2, 95));

hold on ;
smd2(extract_events2==0) = nan;

e2times = zeros(max(extract_events2),1);
for ee = 1: max(extract_events2)
    zvec = zeros(size(smd2));
    zvec(extract_events2 == ee) = smd2(extract_events2 == ee);
    [~, e2times(ee)] = max(zvec);
end
plot(smd2, 'r');
max(extract_events2);

%% prettyplot all subject's time courses together with EB from run 2
figure;
plot([1:450000]./60000,mean(dts1), 'color', [30 83 145]./255, 'LineWidth',1);

hold on;
plot([1:450000]./60000, mean(dts2), 'color', [198 0 27]./255, 'LineWidth',1);
hold on;
[1:450000]./60000;
for ep = 1 : length(e2times)
    line([time_vec(e2times(ep)) time_vec(e2times(ep))],[0 .5 ], 'color', [222 127 139]./255,'LineStyle','-','LineWidth',1);
end


xlim([0 7.333]); ylim([0 0.4]);
xlabel('time (minutes)');
ylabel('agreement');
set(gcf, 'color', 'white');
set(gcf, 'Position', [100 100 1400 200]);

%%

time_vec = [1:450000]./60000;
figure;
subplot(211)
plot([1:450000]./60000,mean(dts1(select_vec,:)), 'color', [131 131 255]./255, 'LineWidth',2);
hold on;
for ep = 1 : length(e1times)
    line([time_vec(e1times(ep)) time_vec(e1times(ep))],[0 .5 ], 'color', 'k','LineStyle','--','LineWidth',1);
end
xlabel('minutes');
ylabel('agreement run 1');

subplot(212)
plot([1:450000]./60000, mean(dts2(select_vec,:)), 'color', [233 53 98]./255, 'LineWidth',2);
hold on;
for ep = 1 : length(e2times)
    line([time_vec(e2times(ep)) time_vec(e2times(ep))],[0 .5 ], 'color', 'k','LineStyle','--','LineWidth',1);
end
xlabel('minutes');
ylabel('agreement run 2');

%
[floor(time_vec(e1times))' round(mod(time_vec(e1times),1)*60)']

[floor(time_vec(e2times))' round(mod(time_vec(e2times),1)*60)']
%


