clear all
close all
clc
path_dataset= '';  % Enter the the path of dataset1 on your system 

path_ = cd;
addpath(genpath([path_ '\function']))
%% Load Dataset

number=1;
for i=1 : 16
    SpikeActStr='SpikeActBZ'+string(i);
    load(SpikeActStr)
    isEmpty=isempty(su);
   
    if(~isEmpty)
        n=size(su,2);
        for j=1 : n
            neurons(number)=su(j);
            number=number+1;
        end    
    end 
    
end
load('Cond')
%% ُRaster Plot
% response after stimulus onset
time_stamps = 1 : 8000;
% response before stimulus onset
time_stamps = time_stamps - 1000;

for i=1: 8
    condition=i;
    ind_h = find(Cond == condition);
    trial_num = 1:size(ind_h, 1);
    figure()
    for j=1:4
        ax = subplot(1, 4, j);
        temp = [];
        temp = neurons{j*j+2}(ind_h, :);
        imagesc(time_stamps, trial_num,temp)           
        colormap(flipud(colormap('gray')));
        xlabel('Time (milisec.)');
        ylabel('Trial number (#)');
        title(strcat('Neuron',num2str(j)));
        text(500, 12, strcat(' N = ', num2str(length(ind_h))),...
       'fontsize', 18, 'fontweight', 'bold', 'Color', 'r');
        text(500, 14, strcat(' Condition = ', num2str(i)),...
       'fontsize', 18, 'fontweight', 'bold', 'Color', 'r');
       line([0 0], ylim, 'Color', 'r');
       ax.XLim = [-1000 7000]; 
       ax.YLim = [0 length(trial_num)];
       set(gca, 'fontsize', 14, 'fontweight', 'bold');
    end    
end


ind_h = find(Cond == 1);
trial_num = 1:size(ind_h, 1);
figure()
ax = subplot(1, 2, 1);
temp = [];
temp = neurons{2}(ind_h, :);
imagesc(time_stamps, trial_num,temp)           
colormap(flipud(colormap('gray')));
xlabel('Time (milisec.)');
ylabel('Trial number (#)');
title(strcat('Condition 1 Angle=1 Radius=1'));
text(500, 12, strcat(' N = ', num2str(length(ind_h))),...
'fontsize', 18, 'fontweight', 'bold', 'Color', 'r');
line([0 0], ylim, 'Color', 'r');
ax.XLim = [-1000 7000]; 
ax.YLim = [0 length(trial_num)];
set(gca, 'fontsize', 14, 'fontweight', 'bold');

ind_h = find(Cond == 9);
trial_num = 1:size(ind_h, 1);
hold on
ax = subplot(1, 2, 2);
temp = [];
temp = neurons{2}(ind_h, :);
imagesc(time_stamps, trial_num,temp)           
colormap(flipud(colormap('gray')));
xlabel('Time (milisec.)');
ylabel('Trial number (#)');
title(strcat('Condition 9 Angle=1 Radius=2'));
text(500, 12, strcat(' N = ', num2str(length(ind_h))),...
'fontsize', 18, 'fontweight', 'bold', 'Color', 'r');
line([0 0], ylim, 'Color', 'r');
ax.XLim = [-1000 7000]; 
ax.YLim = [0 length(trial_num)];
set(gca, 'fontsize', 14, 'fontweight', 'bold');
%% PSTH
win_  = 60;
psth  = @(x) ndass_smooth(1000*mean(x,1), win_);


% response after stimulus onset
t_h = 1 : 8000;
% response before stimulus onset
t_h = t_h - 1000;

randColor = rand(8,3);
for i = 1:4
     figure()
    for j= 1:8
    ind_h = find(Cond == j);
    trial_num = 1 : size(ind_h,1);
    temp = [];
    temp = neurons{i*i+2}(ind_h, :);
    hold on
    plot(t_h, psth(temp), 'Color',randColor(j,:))
    xlabel('Time(milisec.)');
    ylabel('Firing rate (Hz)');
    title(strcat('Neuron', num2str(i)));
    ax.XLim = [-1000 7000];
    ax.YLim = [0 150];
    end
    line([0 0], ylim, 'color', 'r') 
    legend('Condition1','Condition2', 'Condition3', 'Condition4', 'Condition5', 'Condition6', 'Condition7', 'Condition8')

end    
 %% MI (1000-1500)
 % Parameters
win_  = 60;
psth  = @(x) ndass_smooth(1000*mean(x, 1), win_);
T_st  =1;
T_end = 8100;

% Calculation

psth_resp = [];

for tri = 1 : size(Cond, 1)
    for sui = 1 : size(neurons, 2)
        psth_resp(sui, tri, :) = psth(neurons{sui}(tri, :));
    end
end

% Normalization
ind_b = [1 1000];
for  sui = 1 : size(neurons, 2)
    var_h=[squeeze(psth_resp(sui, :, :))];
    var_b=mean2(var_h(:, ind_b(1):ind_b(2)));
    var_max=max(max(var_h));
    fef_psth_resp_norm(sui, :, :)= (squeeze(psth_resp(sui, :, :)) - var_b)/(var_max - var_b);
end

%  Mutual information  signle neuron
win_  = 100;
step_ = 10;
win_h = [1000: step_: 1500; 1000+win_: step_: 1500+win_]';


ind_h1 = find(Cond == 1);
ind_h2 = find(Cond == 2);
ind_h3 = find(Cond== 3);
ind_h4 = find(Cond== 4);
ind_h5 = find(Cond== 5);
ind_h6 = find(Cond== 6);
ind_h7 = find(Cond== 7);
ind_h8 = find(Cond== 8);



total_mi = [];

for sui = 1 : size(neurons,2)
    for ti = 1 : size(win_h,1)
        t1 = win_h(ti,1);
        t2 = win_h(ti,2);
        pref1  = nanmean(fef_psth_resp_norm(sui, ind_h1, t1:t2), 3)';
        pref2  = nanmean(fef_psth_resp_norm(sui, ind_h2, t1:t2), 3)';
        pref3  = nanmean(fef_psth_resp_norm(sui, ind_h3, t1:t2), 3)';
        pref4  = nanmean(fef_psth_resp_norm(sui, ind_h4, t1:t2), 3)';
        pref5  = nanmean(fef_psth_resp_norm(sui, ind_h5, t1:t2), 3)';
        pref6  = nanmean(fef_psth_resp_norm(sui, ind_h6, t1:t2), 3)';
        pref7  = nanmean(fef_psth_resp_norm(sui, ind_h7, t1:t2), 3)';
        pref8  = nanmean(fef_psth_resp_norm(sui, ind_h8, t1:t2), 3)';
       total_mi(sui,ti) = ndass_mi([pref1;pref2;pref3;pref4;pref5;pref6;pref7;pref8],...
        [ones(length(pref1), 1); 2*ones(length(pref2), 1);3*ones(length(pref3),1);4*ones(length(pref4),1);
        5*ones(length(pref5),1);6*ones(length(pref6),1);7*ones(length(pref7),1);8*ones(length(pref8),1)], 20, 5);    
         
    end      
     
    
end

t_h = nanmean(win_h,2)';

figure
hold on
ndass_niceplot(total_mi, t_h, 1, 1, 0, 0)
xlabel('Time (milsec.)');
ylabel('Mutual information  ');
ax.XLim = [1000 1500];
ax.YLim = [0 1];
%line([0 0], ylim, 'color', 'r') 
set(gca, 'fontsize', 14, 'fontweight', 'bold');
%% 2700- 3200
% Raster plot
time_stamps = 2700 : 3200;


for i = 1:4
    figure()
    for j= 1:8
        ind_h = find(Cond == j);
        trial_num = 1 : size(ind_h,1);
        temp = [];
        temp = neurons{i*i+2}(ind_h, :);
        ax = subplot(2, 4, j);
        imagesc(time_stamps, trial_num,temp)           
        colormap(flipud(colormap('gray')));
        xlabel('Time (milisec.)');
        ylabel('Trial number (#)');
        title(strcat('Condition',num2str(j)));
       line([0 0], ylim, 'Color', 'r');
       ax.XLim = [2700 3200]; 
       ax.YLim = [0 length(trial_num)];
       set(gca, 'fontsize', 14, 'fontweight', 'bold');
    end    
end

%PSTH
win_  = 60;
psth  = @(x) ndass_smooth(1000*mean(x,1), win_);


% response after stimulus onset
t_h = 2700 : 3200;
% response before stimulus onset


randColor = rand(8,3);
for i = 1:4
    figure()
    for j= 1:8
        ind_h = find(Cond == j);
        trial_num = 1 : size(ind_h,1);
        temp = [];
        temp = neurons{i*i+2}(ind_h, 2700:3200);
        hold on
        plot(t_h, psth(temp), 'Color',randColor(j,:))
        xlabel('Time(milisec.)');
        ylabel('Firing rate (Hz)');
        title(strcat('Neuron', num2str(i)));
        ax.XLim = [2700 3200];
        ax.YLim = [0 150];
    end
    %line([0 0], ylim, 'color', 'r') 
    legend('Condition1','Condition2', 'Condition3', 'Condition4', 'Condition5', 'Condition6', 'Condition7', 'Condition8')

end    

%MI
win_  = 100;
step_ = 10;
win_h = [2700: step_: 3200; 2700+win_: step_: 3200+win_]';


ind_h1 = find(Cond == 1);
ind_h2 = find(Cond == 2);
ind_h3 = find(Cond== 3);
ind_h4 = find(Cond== 4);
ind_h5 = find(Cond== 5);
ind_h6 = find(Cond== 6);
ind_h7 = find(Cond== 7);
ind_h8 = find(Cond== 8);



t_mi = [];

for sui = 1 : size(neurons,2)
    for ti = 1 : size(win_h,1)
        t1 = win_h(ti,1);
        t2 = win_h(ti,2);
        pref1  = nanmean(fef_psth_resp_norm(sui, ind_h1, t1:t2), 3)';
        pref2  = nanmean(fef_psth_resp_norm(sui, ind_h2, t1:t2), 3)';
        pref3  = nanmean(fef_psth_resp_norm(sui, ind_h3, t1:t2), 3)';
        pref4  = nanmean(fef_psth_resp_norm(sui, ind_h4, t1:t2), 3)';
        pref5  = nanmean(fef_psth_resp_norm(sui, ind_h5, t1:t2), 3)';
        pref6  = nanmean(fef_psth_resp_norm(sui, ind_h6, t1:t2), 3)';
        pref7  = nanmean(fef_psth_resp_norm(sui, ind_h7, t1:t2), 3)';
        pref8  = nanmean(fef_psth_resp_norm(sui, ind_h8, t1:t2), 3)';
       t_mi(sui,ti) = ndass_mi([pref1;pref2;pref3;pref4;pref5;pref6;pref7;pref8],...
        [ones(length(pref1), 1); 2*ones(length(pref2), 1);3*ones(length(pref3),1);4*ones(length(pref4),1);
        5*ones(length(pref5),1);6*ones(length(pref6),1);7*ones(length(pref7),1);8*ones(length(pref8),1)], 20, 5);    
         
    end      
     
    
end

t_h = nanmean(win_h,2)';
%t_h = t_h - 1000;

figure
hold on
ndass_niceplot(t_mi, t_h, 1, 1, 0, 1)
ylabel('Mutual information ');
ax.XLim = [2700 3200];
ax.YLim = [0 1];

set(gca, 'fontsize', 14, 'fontweight', 'bold');



FirstMIMean=nanmean(total_mi,2);
SecondMIMean=nanmean(t_mi,2);
scatter(FirstMIMean,SecondMIMean,'filled');
xlabel('MI in 1000-1500 ms');
ylabel('MI in 2700-3200 ms');
[h,p]=ttest(FirstMIMean,SecondMIMean);


 %% ROC
win_  = 100;
step_ = 10;
%win_h = [1000: step_: 1500; win_+1000: step_: 1500+win_]'; % for [1000-1500]
win_h = [2700: step_: 3200; win_+2700: step_: 3200+win_]'; % for[2700-3200]
ind_h1 = find(Cond== 1);
ind_h2 = find(Cond== 2);
ind_h9 = find(Cond== 9);
ind_h10 = find(Cond== 10);

 
for sui = 1 : size(neurons,2)
    for ti = 1 : size(win_h,1)
        t1 = win_h(ti,1);
        t2 = win_h(ti,2);
        pref1  = nanmean(fef_psth_resp_norm(sui, ind_h1, t1:t2), 3)';
        pref2  = nanmean(fef_psth_resp_norm(sui, ind_h9, t1:t2), 3)';
        
       roc1(sui,ti) = ndass_roc( pref1,pref2);
         
        pref1  = nanmean(fef_psth_resp_norm(sui, ind_h2, t1:t2), 3)';
        pref2  = nanmean(fef_psth_resp_norm(sui, ind_h10, t1:t2), 3)';
        
       roc2(sui,ti) = ndass_roc( pref1,pref2);
    end      
     
    
end

t_h = nanmean(win_h, 2)';
figure
hold on
ndass_niceplot(roc1, t_h, 1, 0, 1, 1)
ndass_niceplot(roc2, t_h, 1, 0, 0, 1)
xlabel('Time(milisec.)');
ylabel('AROC (a.u.)');
ax.XLim = [1000 1500];
ax.YLim = [0 1];
set(gca, 'fontsize', 14, 'fontweight', 'bold');

%% SVM
win_  = 100;
step_ = 10;
win_h = [1000: step_: 1500; win_+1000: step_: 1500+win_]';


ind_h1 = find(Cond == 1);
ind_h2 = find(Cond == 2);
ind_h3 = find(Cond== 3);
ind_h4 = find(Cond== 4);
ind_h5 = find(Cond== 5);
ind_h6 = find(Cond== 6);
ind_h7 = find(Cond== 7);
ind_h8 = find(Cond== 8);


preformance = [];

for ti = 1 : size(win_h, 1)
    t1 = win_h(ti, 1);
    t2 = win_h(ti, 2);
    
    group1 = nanmean(fef_psth_resp_norm(:, ind_h1, t1:t2), 3)';
    group2 = nanmean(fef_psth_resp_norm(:, ind_h2, t1:t2), 3)';
    group3 = nanmean(fef_psth_resp_norm(:, ind_h3, t1:t2), 3)';
    group4 = nanmean(fef_psth_resp_norm(:, ind_h4, t1:t2), 3)';
    group5 = nanmean(fef_psth_resp_norm(:, ind_h5, t1:t2), 3)';
    group6 = nanmean(fef_psth_resp_norm(:, ind_h6, t1:t2), 3)';
    group7= nanmean(fef_psth_resp_norm(:,  ind_h7, t1:t2), 3)';
    group8 = nanmean(fef_psth_resp_norm(:, ind_h8, t1:t2), 3)';

    svm_res = ndass_svm([group1; group2; group3;group4;group5;group6;group7;group8],...
        [ones(length(group1),1); 2*ones(length(group2),1); 3*ones(length(group3),1);4*ones(length(group4),1);5*ones(length(group5),1);6*ones(length(group6),1);7*ones(length(group7),1);8*ones(length(group8),1)], ...
        0.7, 10) ;
    preformance(:,ti) = ((svm_res.pt)-1/8)/(1-1/8);
    
end

t_h = nanmean(win_h, 2)';

figure
hold on
ndass_niceplot(preformance, t_h, 1, 1, 0, 0)
xlabel('Time(milisec.)');
ylabel('Performance ');
ax.XLim = [1000 1500];
ax.YLim = [0 1];
set(gca, 'fontsize', 14, 'fontweight', 'bold');
%% LFP %%
%% Load Dataset
load('LFP16')
lfp1=LFP;
load('LFP1')
lfp2=LFP;
load('Cond');

% Remove Artifact
artif1=ndass_rmartifact(lfp1,1,3);
lfp1 = ndass_nothfilter(lfp1,50);
nomalizer = @(x) (x - nanmean(x(:)))/(nanstd(x(:)));
lfp1 = nomalizer(lfp1);

artif2=ndass_rmartifact(lfp2,1,3);
lfp2 = ndass_nothfilter(lfp2,50);
nomalizer = @(x) (x - nanmean(x(:)))/(nanstd(x(:)));
lfp2 = nomalizer(lfp2);
t1 = 1000; t2 = 8000;
tb1 = 1; tb2 = 1000;
%% Time Power
% multitaper
chronux=true;
params.Fs=1000; % sampling frequency
params.fpass=[1 100]; % band of frequencies to be kept
params.tapers=[3 5];%[3 5]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[0 0.05];
params.trialave=0;
movingwindow=[0.2 0.05];

[S,t,f]=mtspecgramc(lfp1',movingwindow,params);
figure
for i=1:8
    ind_h = ismember(Cond,i)&(~artif1);
    subplot(2,4,i)
    temp=[];
    temp=squeeze(nanmean(S(:,:,ind_h),3)');
    h=pcolor(t-1,f,temp);
    h.EdgeColor='none';
    colormap(jet);
    clb=colorbar;
    title(strcat('Condition',num2str(i)));
    xlabel('time (sec)');
    ylabel('power');
end



% wavelet
fs = 1000;
freqs2use = [5:40 40:10:100];
[analytic_sig, freqs2use] = ndass_wavelet(lfp1, freqs2use, fs);

t =([1:8000])/1000;
S = abs(analytic_sig).^2;
figure()
for i=1:8
    ind_h = ismember(Cond,i)&(~artif1);
    subplot(2,4,i)
    temp=[];
    temp=squeeze(nanmean(S(ind_h,:,:),1));
    h=pcolor(t-1,freqs2use,temp);
    h.EdgeColor='none';
    colormap(jet);
    clb=colorbar;
    title(strcat('Condition',num2str(i)));
    xlabel('time (sec)');
    ylabel('power');
end
%% Phase Phase locking
% over trial
artif  = (artif1 | artif2);
ind_h1 = ismember(Cond,1)&(~artif);
ind_h2 = ismember(Cond,2)&(~artif);
ind_h3 = ismember(Cond,3)&(~artif);
ind_h4 = ismember(Cond,4)&(~artif);
ind_h5 = ismember(Cond,5)&(~artif);
ind_h6 = ismember(Cond,6)&(~artif);
ind_h7 = ismember(Cond,7)&(~artif);
ind_h8 = ismember(Cond,8)&(~artif);

fs = 1000;
t  = [1:8000]/1000; % in second
f  = [5:40 40:10:100];
disp ('wavelet calculation take time')
[analytic_sig1, f] = ndass_wavelet(lfp1, f, fs);
[analytic_sig2, f] = ndass_wavelet(lfp2, f, fs);
disp ('finish wavelet calculation ')

ppl = @ (x1,x2) abs(nanmean(exp(1i*(x1-x2))));

phi1 = angle(analytic_sig1);
phi2 = angle(analytic_sig2);
figure( 'Name','PPL LFP1&LFP16')
for i=1:8
    ind_h = ismember(Cond,i)&(~artif);
    subplot(2,4,i)
    temp = squeeze(ppl(phi1(ind_h,:,:),phi2(ind_h,:,:)));
    var_b = repmat(nanmean(temp(:,t>0&t<0.1000),2),1,length(t));
    std_base = repmat(nanstd(temp(:,t>0&t<0.1000),[],2),1,length(t));
    temp = (temp -var_b)./std_base ;
    h=pcolor(t-1,f,temp);
    h.EdgeColor = 'none';
    colormap(jet); 
    clb= colorbar;
    title (strcat('PPL for condition',num2str(i)));
    xlabel('Time from sample onset ');
    ylabel('Frequency (Hz)');
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
end
%% PPL Over time (1000-1500)
phi1 = angle(analytic_sig1);
phi2 = angle(analytic_sig2);
trialnumber=size(phi1,1);  freqnumber=size(phi1,2); timepointnumber=size(phi1,3);
figure( 'Name','PPL')%figure('Position',[187  558  1676  420], 'Name','PLV')
timebin=[2700:3200];
timebinbase=[1:1000];
figure( 'Name','PPL LFP1&LFP16')
for i=1:8
    ind_h = ismember(Cond,i)&(~artif);
    subplot(2,4,i)
    trialnums=find(ind_h==1);
    var_h=[];var_b1=[];
    for fi =1:freqnumber
    for tr = 1 :length(trialnums) 
         var_h(fi,tr)=ndass_ppl(phi1(tr,fi,timebin),phi2(tr,fi,timebin));
         var_b1(fi,tr)=ndass_ppl(phi1(tr,fi,timebinbase),phi2(tr,fi,timebinbase));      
    end
    end
    plot(f,nanmean(var_h,2),'color',[0.5,1,0]); hold on; plot(f,nanmean(var_b1,2),'color',[0.1,0.5,1]); xlabel('frequency'); ylabel('PPL');
    ndass_niceplot(var_h',f, 1, 0.5, 1, 1); ndass_niceplot(var_b1',f, 1, 0.1, 0.5, 1); 
    var_b = repmat(nanmean(var_b1,2),1,size(var_b1,2));
    std_base = repmat(nanstd(var_b1,[],2),1,size(var_b1,2));
    var_h = (var_h -var_b)./std_base ;
    ndass_niceplot(var_h',f, 1, 0.5, 0.5, 0); 
    title (strcat('PPL for condition ',num2str(i)));
    ylabel('PPL ');
    xlabel('Frequency (Hz)');
    set(gca, 'fontsize', 14, 'fontweight', 'bold');
end

%% Coherency multitapper

t1 = 1000; t2 = 1500; 
% Settings Chronux
params.Fs=1000; % sampling frequency
params.fpass=[1 100]; % band of frequencies to be kept
params.tapers=[3 5];%[3 5]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[0 0.05];
params.trialave=1;
S1_mt = []; S2_mt = [];
[C1,phi1,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h1,t1:t2)',lfp2(ind_h1,t1:t2)',params);  % Attention: 't' is in seconds!
[C2,phi2,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h2,t1:t2)',lfp2(ind_h2,t1:t2)',params);  % Attention: 't' is in seconds!
[C3,phi3,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h3,t1:t2)',lfp2(ind_h3,t1:t2)',params);  % Attention: 't' is in seconds!
[C4,phi4,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h4,t1:t2)',lfp2(ind_h4,t1:t2)',params);  % Attention: 't' is in seconds!
[C5,phi5,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h5,t1:t2)',lfp2(ind_h5,t1:t2)',params);  % Attention: 't' is in seconds!
[C6,phi6,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h6,t1:t2)',lfp2(ind_h6,t1:t2)',params);  % Attention: 't' is in seconds!
[C7,phi7,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h7,t1:t2)',lfp2(ind_h7,t1:t2)',params);  % Attention: 't' is in seconds!
[C8,phi8,S12,S1,S2,f_mt] = coherencyc(lfp1(ind_h8,t1:t2)',lfp2(ind_h8,t1:t2)',params);  % Attention: 't' is in seconds!

randColor = rand(8,3);
figure( 'Name','Coherency')
for i=1:8
     name=strcat('C',num2str(i));
     var_h = []; 
     var_h = eval(name)';
     hold on
     plot(f_mt,var_h,'Color',randColor(i,:))
     title ('Coherency')
     xlabel('Frequency (Hz)');
     ylabel('Coherency (a.u.)');
     set(gca, 'fontsize', 14, 'fontweight', 'bold');
end
legend('C1','C2','C3','C4','C5','C6','C7','C8');

% Time Ferquence
movingwin=[0.2 0.050];%movingwin=[0.5 0.05];

[C1,phi1,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h1,:)',lfp2(ind_h1,:)',movingwin,params);  % Attention: 't' is in seconds!
[C2,phi2,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h2,:)',lfp2(ind_h2,:)',movingwin,params);  % Attention: 't' is in seconds!
[C3,phi3,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h3,:)',lfp2(ind_h3,:)',movingwin,params);  % Attention: 't' is in seconds!
[C4,phi4,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h4,:)',lfp2(ind_h4,:)',movingwin,params);  % Attention: 't' is in seconds!
[C5,phi5,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h5,:)',lfp2(ind_h5,:)',movingwin,params);  % Attention: 't' is in seconds!
[C6,phi6,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h6,:)',lfp2(ind_h6,:)',movingwin,params);  % Attention: 't' is in seconds!
[C7,phi7,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h7,:)',lfp2(ind_h7,:)',movingwin,params);  % Attention: 't' is in seconds!
[C8,phi8,S12,S1,S2,t_mt,f_mt] = cohgramc(lfp1(ind_h8,:)',lfp2(ind_h8,:)',movingwin,params);  % Attention: 't' is in seconds!


figure( 'Name','Coherency')
for i=1:8
     subplot(2,4,i)
     name=strcat('C',num2str(i));
     var_h = []; 
     var_h = eval(name)';
     var_b = repmat(nanmean(var_h(:,t_mt>0&t_mt<1),2),1,length(t_mt));
     std_base = repmat(nanstd(var_h(:,t_mt>0&t_mt<1),[],2),1,length(t_mt));
     var_h = (var_h -var_b)./std_base ;
     h=pcolor(t_mt-1,f_mt,var_h);
     h.EdgeColor = 'none';
     colormap(jet); 
     clb= colorbar;
     title (strcat('Coherency for condition',num2str(i)))
     xlabel('Time (sec) ');
     ylabel('Frequency (Hz)');
     set(gca, 'fontsize', 14, 'fontweight', 'bold');
end
%% Spike Field Coherency
t1 = 1000; t2 = 1500;
t = 1: 8000;
ind_t = t1<t&t<t2;

number=1;
for i=1 : 16
    load(strcat('LFP',num2str(i)));
    lfps{i}=LFP; 
end   

%First neuron spike
Sample_neuron=neurons{2};
lfp1=lfps{4};
lfp2=lfps{5};
lfp3=lfps{6};
artif1 = ndass_rmartifact(lfp1,0,3);
artif2 = ndass_rmartifact(lfp2,0,3);
artif3 = ndass_rmartifact(lfp3,0,3);
artif  = (artif1 | artif2 | artif3);

% Notch filter
lfp1 = ndass_nothfilter(lfp1,50);
lfp2 = ndass_nothfilter(lfp2,50);
lfp3 = ndass_nothfilter(lfp3,50);

% Normalizaton
nomalizer = @(x) (x - mean2(x))/nanstd(x(:));
 lfp1 = nomalizer(lfp1);
 lfp2 = nomalizer(lfp2);
 lfp3 = nomalizer(lfp3);

ind_h1 = ismember(Cond,1)&(~artif);
ind_h2 = ismember(Cond,2)&(~artif);
ind_h3 = ismember(Cond,3)&(~artif);
ind_h4 = ismember(Cond,4)&(~artif);
ind_h5 = ismember(Cond,5)&(~artif);
ind_h6 = ismember(Cond,6)&(~artif);
ind_h7 = ismember(Cond,7)&(~artif);
ind_h8 = ismember(Cond,8)&(~artif);

params.Fs=1000; % sampling frequency
params.fpass=[1 100]; % band of frequencies to be kept
params.tapers=[3 5];%[3 5]; % taper parameters
params.pad=1; % pad factor for fft
params.err=[1 0.05];
params.trialave = 1;
randColor = rand(8,3);
figure( 'Name','Spike Field Coherency with Neuron 2 Spike')
subplot(131);
for i=1:8
    hold on
    ind_h=find(Cond==i);
    datalfp = lfp1(ind_h,ind_t)';
    sp = Sample_neuron(ind_h,ind_t);
    datasp =[];
    for j = 1:size(sp,1)
         datasp(j).sp = find(sp(j,:))/1000;
    end

   [C,phi,S12,S1,S2,f_mt,zerosp,confC,phistd]=coherencycpt(datalfp,datasp,params);
   plotsig(C,confC,1,f_mt,'Color',randColor(i,:)); 
   set(gca,'FontName','Times New Roman','Fontsize', 16);
   xlabel(''); ylabel('Coherence');
end

subplot(132);
for i=1:8
    hold on
    ind_h=find(Cond==i);
    datalfp = lfp2(ind_h,ind_t)';
    sp = Sample_neuron(ind_h,ind_t);
    datasp =[];
    for j = 1:size(sp,1)
         datasp(j).sp = find(sp(j,:))/1000;
    end

   [C,phi,S12,S1,S2,f_mt,zerosp,confC,phistd]=coherencycpt(datalfp,datasp,params);
   plotsig(C,confC,1,f_mt,'Color',randColor(i,:)); 
   set(gca,'FontName','Times New Roman','Fontsize', 16);
   xlabel(''); ylabel('Coherence');
end
subplot(133);
for i=1:8
    hold on
    ind_h=find(Cond==i);
    datalfp = lfp3(ind_h,ind_t)';
    sp = Sample_neuron(ind_h,ind_t);
    datasp =[];
    for j = 1:size(sp,1)
         datasp(j).sp = find(sp(j,:))/1000;
    end

   [C,phi,S12,S1,S2,f_mt,zerosp,confC,phistd]=coherencycpt(datalfp,datasp,params);
   plotsig(C,confC,1,f_mt,'Color',randColor(i,:)); 
   set(gca,'FontName','Times New Roman','Fontsize', 16);
   xlabel(''); ylabel('Coherence');
end