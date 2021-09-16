close all
clear variables
clc

%% Import data from TDT files
% Update with current folder path
data = TDTbin2mat(uigetdir);

%% Create time vector from data
% Fs = sampling rate
% tssub = take out first 10 seconds of data

ts = (0:(length(data.streams.x65B.data)-1))/data.streams.x65B.fs;
Fs=data.streams.x65B.fs;
tssub=ts(round(10*Fs):end);

%% Extract signal and control data and create plot
% Apply lowpass filter 

sig=data.streams.x65B.data(round(10*Fs):end);
ctr=data.streams.x05V.data(round(10*Fs):end);
sig=double(sig);
ctr=double(ctr);
ctr_filt=lowpassphotometry(ctr,Fs,2);

plot(tssub,sig)
hold on
plot(tssub,ctr_filt)
xlabel('Time(s)','FontSize',14)
ylabel('mV at detector','FontSize',14)
daspect([1 1 1])
hold off
saveas(gcf, 'Raw data figure.emf')

%% Normalize signal to control channel, plot and save figure
% Add notes to plot

[normDat] = deltaFF(sig,ctr_filt);
figure;plot(tssub, normDat)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('deltaF/F','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')

saveas(gcf, 'Figure dFF and Notes.emf')

hold off
