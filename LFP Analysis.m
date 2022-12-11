clear all
close all
clc
path_ = cd;
addpath(genpath([path_ '\function']))

%% Load Dataset
load('LFP16')
lfp1=LFP;
load('LFP1')
lfp2=LFP;
load('Cond');
%%
artif1=ndass_rmartifact(lfp1,1,3);
lfp1 = ndass_nothfilter(lfp1,50);
nomalizer = @(x) (x - nanmean(x(:)))/(nanstd(x(:)));
lfp1 = nomalizer(lfp1);

artif2=ndass_rmartifact(lfp2,1,3);
lfp2 = ndass_nothfilter(lfp2,50);
nomalizer = @(x) (x - nanmean(x(:)))/(nanstd(x(:)));
lfp2 = nomalizer(lfp2);

ind_h1 = ismember(Cond,1)&(~artif1);
ind_h2 = ismember(Cond,2)&(~artif1);
ind_h3 = ismember(Cond,3)&(~artif1);
ind_h4 = ismember(Cond,4)&(~artif1);
ind_h5 = ismember(Cond,5)&(~artif1);
ind_h6 = ismember(Cond,6)&(~artif1);
ind_h7 = ismember(Cond,7)&(~artif1);
ind_h8 = ismember(Cond,8)&(~artif1);

t1 = 1000; t2 = 8000;
tb1 = 1; tb2 = 1000;
% multitaper
chronux=true;
params.Fs=1000; % sampling frequency
params.fpass=[1 100]; % band of frequencies to be kept
params.tapers=[2 3];%[3 5]; % taper parameters
params.pad=2; % pad factor for fft
params.err=[0 0.05];
params.trialave=0;
movingwindow=[0.2 0.05];

[S,t,f]=mtspecgramc(lfp1',movingwindow,params);
%[S1,f]=mtspectrumc(lfp1(ind_h1,t1:t2)',params);
%S1 = S1';
%[S4,f]=mtspectrumc(lfp1(ind_h4,tb1:tb2)',params);
%S4= S4';

figure
subplot(2,4,1)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h1),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('condition1');
xlabel('time');
ylabel('power');

subplot(2,4,2)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h2),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 2');
xlabel('time');
ylabel('power');

subplot(2,4,3)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h3),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 3');
xlabel('time');
ylabel('power');

subplot(2,4,4)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h4),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 4');
xlabel('time');
ylabel('power');

subplot(2,4,5)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h5),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 5');
xlabel('time');
ylabel('power');

subplot(2,4,6)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h6),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 6');
xlabel('time');
ylabel('power');

subplot(2,4,7)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h7),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 7');
xlabel('time');
ylabel('power');

subplot(2,4,8)
temp=[];
temp=squeeze(nanmean(S(:,:,ind_h8),3)');
h=pcolor(t-0.5,f,temp);
h.EdgeColor='none';
colormap(jet);
clb=colorbar;
title('Condition 8');
xlabel('time');
ylabel('power');




%%
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


trialnumber=size(phi1,1);  freqnumber=size(phi1,2); timepointnumber=size(phi1,3);
figure( 'Name','PLV')
subplot(131)
var_h = []; 

var_h = squeeze(ppl(phi1(ind_h1,:,:),phi2(ind_h1,:,:)));

var_b = repmat(nanmean(var_h(:,t>0&t<0.1000),2),1,length(t));
std_base = repmat(nanstd(var_h(:,t>0&t<0.1000),[],2),1,length(t));
var_h = (var_h -var_b)./std_base ;
h=pcolor(t-1,f,var_h);
h.EdgeColor = 'none';
colormap(jet); 
clb= colorbar;
title ('PLV for condition1')
xlabel('Time from sample onset ');
ylabel('Frequency (Hz)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
caxis([-3 3])
% % % a sample polar plot
% % % time=2000; fr=10; phasediff=squeeze(phi1(ind_h1,fr,time)-phi2(ind_h1,fr,time));
% % % phasediff=phasediff'; figure; polar([zeros(size(phasediff));phasediff],[zeros(size(phasediff)); ones(size(phasediff))]);

% a faster way
var_h = squeeze(ppl(phi1(ind_h1,:,:),phi2(ind_h1,:,:)));

var_b = repmat(nanmean(var_h(:,t>0&t<0.500),2),1,length(t));
std_base = repmat(nanstd(var_h(:,t>0&t<0.500),[],2),1,length(t));
var_h = (var_h -var_b)./std_base ;
h=pcolor(t-0.5,f,var_h);
h.EdgeColor = 'none';
colormap(jet); 
clb= colorbar;
title ('PLV for condition1')
xlabel('Time from sample onset ');
ylabel('Frequency (Hz)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
caxis([-3 3])

subplot(132)
var_h = []; 
var_h = squeeze(ppl(phi1(ind_h2,:,:),phi2(ind_h2,:,:)));
var_b = repmat(nanmean(var_h(:,t>0&t<0.500),2),1,length(t));
std_base = repmat(nanstd(var_h(:,t>0&t<0.500),[],2),1,length(t));
var_h = (var_h -var_b)./std_base ;
h=pcolor(t-0.5,f,var_h);
h.EdgeColor = 'none';
colormap(jet); 
clb= colorbar;
title ('PLV for condition2')
xlabel('Time from sample onset ');
ylabel('Frequency (Hz)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
caxis([-3 3])

subplot(133)
var_h = []; 
var_h = squeeze(ppl(phi1(ind_h3,:,:),phi2(ind_h3,:,:)));
var_b = repmat(nanmean(var_h(:,t>0&t<0.500),2),1,length(t));
std_base = repmat(nanstd(var_h(:,t>0&t<0.500),[],2),1,length(t));
var_h = (var_h -var_b)./std_base ;
h=pcolor(t-0.5,f,var_h);
h.EdgeColor = 'none';
colormap(jet); 
clb= colorbar;
title ('PLV for condition 3')
xlabel('Time from sample onset ');
ylabel('Frequency (Hz)');
set(gca, 'fontsize', 14, 'fontweight', 'bold');
caxis([-3 3])