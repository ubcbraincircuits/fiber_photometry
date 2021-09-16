close all
clear variables
clc

%% Import data from TDT files
% Update with current folder path
data = TDTbin2mat(uigetdir);
%filename = [ 'July-' data.info.blockname(1,:)]
filename = [data.info.blockname(1,:)];


%% Create time vector from data
% Fs = sampling rate
% tssub = take out first 10 seconds of data

ts = (0:(length(data.streams.x65B.data)-1))/data.streams.x65B.fs;
Fs=data.streams.x65B.fs;
TotalTime = length(ts)/Fs;

if TotalTime < 570
    disp('Warning: Total Time is less than 10 minutes')
else
    tssub = ts(round(30*Fs):round(600*Fs));
end

%Check length of tssub
lengthtime = length(tssub)/Fs;
%% Extract signal and control data and create plot
% Apply lowpass filter 

sig=data.streams.x65B.data(round(30*Fs):round(600*Fs));
ctr=data.streams.x05V.data(round(30*Fs):round(600*Fs));
sig = double(sig);
ctr = double(ctr);

ctr_filt=lowpassphotometry(ctr,Fs,2);

plot(tssub,sig)
hold on
plot(tssub,ctr_filt-100)
xlabel('Time(s)','FontSize',14)
ylabel('mV at detector','FontSize',14)
ylim([0 700]);
hold off
saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Raw Data Graphs\',filename,'.emf']);

%% Normalize signal to control channel, plot and save figure
% Add notes to plot

[normDat] = deltaFF(sig,ctr_filt);
%normDatfill = filloutliers(normDat, 'linear');
%normDat(round(150*Fs):round(250*Fs)) = normDatfill(round(150*Fs):round(250*Fs));
figure;plot(tssub, normDat)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('deltaF/F','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
ylim([-10 25])
hold off
saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\DFF Graphs\',filename,'.emf']);

% %% Calculate mean, SD, and z score
% Mean = mean(normDat);
% SD = std(normDat);
% zscore = normDat/SD;
% figure;plot(tssub, normDat)
% daspect([1 0.1 1]);
% hold on
% xlabel('Time(s)','FontSize',14)
% ylabel('Z Score','FontSize',14)
% title(data.info.blockname(1,:),'Interpreter','none')
% hold off
% saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Z score Graphs\',filename,'.emf']);
% 
% zscoremean = mean(zscore);
% disp(Mean);
% disp(SD);
% disp(zscoremean);

% %% Find peaks with unsmoothed data, plot and save
% % Smooth z score data first
% % Peak prominence threshold is 1 STD of normDat
% 
% figure; findpeaks(zscore, 'MinPeakProminence', 0.5, 'annotate','extents');
% totalpeaks = numel(findpeaks(zscore, 'MinPeakProminence', 0.5, 'annotate','extents'));
% [pks,locs,w,p] = findpeaks(zscore, 'MinPeakProminence', 0.5,'annotate','extents');
% allpks=[locs;pks;w;p]';
% peaklocs_time=(locs);
% peaklocs_time=peaklocs_time';
% hold on
% ylabel('Z Score','FontSize',14)
% title(data.info.blockname(1,:),'Interpreter','none')
% daspect([1 0.0001 1])
% saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Peak Plots\',filename,'.emf']);
% hold off 
% 
% save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\All Peaks Info\', filename],'allpks', '-ascii');
% %% Calculate frequency and save 
% % Adjust peak prominence as needed
% 
% PeakCount = totalpeaks;
% PeakFreq = PeakCount/lengthtime;
% PeakFreqpermin = PeakFreq * 60;
% 
% disp('Peak Freq per min is');
% disp(PeakFreqpermin);
% 
% save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Count\', filename],'PeakCount', '-ascii');
% save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Freq per min\', filename],'PeakFreqpermin', '-ascii');
% 
% 
% %% Calculate average peak amplitude and width for individual epocs and save
% 
% figure; plot(peaklocs_time, p, '-o', 'MarkerFaceColor', 'b');
% hold on
% xlabel('Time(s)','FontSize',14)
% ylabel('Peak amplitude (ZdF/F)','FontSize',14)
% hold off
% saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Peak Amp Plots\',filename,'.emf']);
% 
% width = (w/Fs) *1000;
% figure; plot(peaklocs_time, width, '-o', 'MarkerFaceColor', 'b');
% hold on
% xlabel('Time(s)','FontSize',14)
% ylabel('Peak width (s)','FontSize',14)
% hold off
% saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Peak Width Plots\',filename,'.emf']);
% 
% PeakAmpAverage = mean(p);
% PeakWidthAverage = mean(width);
% disp('Average Peak Amp is');
% disp(PeakAmpAverage);
% disp('Average Peak Width is');
% disp(PeakWidthAverage);
% 
% save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Amplitudes\', filename],'PeakAmpAverage', '-ascii');
% save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Widths\', filename],'PeakWidthAverage', '-ascii');
% 

%% Calculate frequency of smoothed data and save 
% Adjust peak prominence as needed

smoothnormDat = smooth(normDat,Fs);
smoothnormDat = smoothnormDat';
smoothzscore = smoothnormDat/(std(smoothnormDat));
filename = [filename '-smoothed'];

figure; findpeaks(smoothzscore, 'MinPeakProminence', 0.5, 'annotate','extents');
totalpeakssmooth = numel(findpeaks(smoothzscore, 'MinPeakProminence', 0.5, 'annotate','extents'));
[pks,locs,w,p] = findpeaks(smoothzscore, 'MinPeakProminence', 0.5,'annotate','extents');
allpks2=[locs;pks;w;p]';
peaklocs_time2=(locs);
peaklocs_time2=peaklocs_time2';
hold on
ylabel('Z Score','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
daspect([1 0.00005 1])
saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Peak Plots (smoothed)\',filename,'.emf']);
hold off 

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\All Peaks Info (smoothed)\', filename],'allpks2', '-ascii');

PeakCount2 = totalpeakssmooth;
PeakFreq2 = totalpeakssmooth/lengthtime;
PeakFreqpermin2 = PeakFreq2 * 60;

disp('Peak Freq per min smoothed is');
disp(PeakFreqpermin2);

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Count (smoothed)\', filename],'PeakCount2', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Freq per min (smoothed)\', filename],'PeakFreqpermin2', '-ascii');


%% Calculate average peak amplitude and width for individual epocs and save

figure; plot(peaklocs_time2, p, '-o', 'MarkerFaceColor', 'b');
hold on
xlabel('Time(s)','FontSize',14)
ylabel('Peak amplitude (mv)','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
hold off
saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Peak Amp Plots (smoothed)\',filename,'.emf']);


width2 = (w/Fs) *1000;
figure; plot(peaklocs_time2, w, '-o', 'MarkerFaceColor', 'b');
hold on
xlabel('Time(s)','FontSize',14)
ylabel('Event width (ms)','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
hold off
saveas(gcf,['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Peak Width Plots (smoothed)\',filename,'.emf']);


PeakAmpAverage2 = mean(p);
PeakWidthAverage2 = mean(width2);
disp('Average Peak Amp is');
disp(PeakAmpAverage2);
disp('Average Peak Width is');
disp(PeakWidthAverage2);

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Amplitudes\', filename],'PeakAmpAverage2', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\Open Field\Open Field - 6-7 month GCAMP cohort\Overall Peak Widths\', filename],'PeakWidthAverage2', '-ascii');


%%
%DatatoCopy = [Mean; SD; zscoremean; PeakFreqpermin; PeakAmpAverage; PeakWidthAverage; PeakFreqpermin2; PeakAmpAverage2; PeakWidthAverage2];
