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

sig=data.streams.x65G.data(round(10*Fs):end);
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
%saveas(gcf, 'Raw data figure.emf')

%% Normalize signal to control channel, plot and save figure
% Add notes to plot

[normDat] = deltaFF(sig,ctr_filt);
figure;plot(tssub, normDat)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('deltaF/F','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
vline([data.epocs.Note.onset(1,1) data.epocs.Note.onset(2,1) data.epocs.PrtA.onset(1,1) data.epocs.PrtA.onset(2,1)], {'k', 'k', 'g', 'g'});

%saveas(gcf, 'Figure dFF and Notes.emf')

hold off

%%
Notepoints = (data.epocs.Note.onset.*Fs);
Notepoints = Notepoints - round(10*Fs);
Rotarodpoints = (data.epocs.PrtA.onset.*Fs);
Rotarodpoints = Rotarodpoints - round(10*Fs); 

%Notesoff = (data.epocs.Note.offset.*Fs);
%Notesoff = Notesoff - 10000;
Pickup = round(Notepoints(1,1));
Start = round(Rotarodpoints(1,1));
Stop = round(Rotarodpoints(2,1));
Down = round(Notepoints(2,1));
%End = round(Notesoff(4,1));

%% Calculate New z score using baseline and create plot
Before = normDat((Pickup-round(280*Fs):(Pickup-round(10*Fs))));
Up = normDat(Pickup:Start);
Rotarod = normDat(Start:Stop);
Putdown = normDat(Stop:Down);
After = normDat(Down+round(10*Fs):Down+round(280*Fs));

stdev=std(Before);
zscore=normDat/stdev;
%tssub = tssub - data.epocs.PrtA.onset(1,1);

figure;plot(tssub(300*Fs:600*Fs), zscore(300*Fs:600*Fs));
hold on
xlabel('Time(s)','FontSize',14)
ylabel('Z Score','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
vline([data.epocs.Note.onset(1,1) data.epocs.Note.onset(2,1) data.epocs.PrtA.onset(1,1) data.epocs.PrtA.onset(2,1)], {'k', 'k', 'g', 'g'});

saveas(gcf, 'New Z Score and Notes.emf')

hold off

%% Downsample data

filt_zscore=lowpassphotometry(zscore,Fs,1);
zscore_down=downsample(filt_zscore,1000)';
ts_down=downsample(tssub,1000)';
%figure; plot(ts_down, zscore_down)

%% Set "start" as 0, plot, save figure and raw data
% Second note should be "start" 
tssub_start = tssub - data.epocs.PrtA.onset(1,1);
note_start = data.epocs.Note.onset - data.epocs.PrtA.onset(1,1);
rotarod_stop = data.epocs.PrtA.onset(2,1) - data.epocs.PrtA.onset(1,1);

figure;plot(tssub_start, zscore)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('Z Score','FontSize',14)
title(data.info.blockname(1,:),'Interpreter','none')
vline([note_start(1,1) note_start(2,1) 0 rotarod_stop(1,1)], {'k', 'k', 'g', 'g'});
daspect([1 0.05 1])
hold off
saveas(gcf, 'Figure Z Score Start Align.emf');

tssub_start_down=downsample(tssub_start,1000)';

startalign = NaN(300,3);
startalign(1:length(tssub_start_down),1) = tssub_start_down;
startalign(1:length(zscore_down),2) = zscore_down;
startalign(1:length(note_start),3) = note_start;
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New Start Align\', data.info.blockname(1,:)],'startalign', '-ascii');

%% Calculate and save new Z score means 

ZBefore = zscore(Pickup-round(280*Fs):(Pickup-round(10*Fs)));
ZPickup = zscore(Pickup:Start);
ZRotarod = zscore(Start:Stop);
ZPutdown = zscore(Stop:Down);
ZAfter = zscore(Down+round(10*Fs):Down+round(280*Fs));

Beforemean = mean(ZBefore);
Pickupmean = mean(ZPickup);
Trialmean = mean(ZRotarod);
Putdownmean = mean(ZPutdown);
Aftermean = mean(ZAfter);

Meansnew = [Beforemean Pickupmean Trialmean Putdownmean Aftermean];
figure; plot(Meansnew)
Meanshift = Meansnew - Beforemean
figure; plot(Meanshift)

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New Z Scores\', data.info.blockname(1,:)],'Meansnew', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New Norm Z Scores\', data.info.blockname(1,:)],'Meanshift', '-ascii');

%% Calculate mean for 30 and 60s intervals up to 2 minutes on rotarod

if Start+round(30*Fs) > Stop
    thirty1 = NaN
else
thirty1 = mean(zscore(Start:Start+round(30*Fs)));
end

if Start+round(60*Fs) > Stop
    thirty2 = NaN
else
thirty2 = mean(zscore(Start+(round(30*Fs)+1):Start+round(60*Fs)));
end

if Start + round(90*Fs) > Stop
    thirty3 = NaN
else
thirty3 = mean(zscore(Start+(round(60*Fs)+1):Start+round(90*Fs)));
end

if Start + round(120*Fs) > Stop
    thirty4 = NaN
else
thirty4 = mean(zscore(Start+(round(90*Fs)+1):Start+round(120*Fs)));
end

thirty = [thirty1 thirty2 thirty3 thirty4];
figure; plot(thirty)

if Start + round(60*Fs) > Stop
    sixty = NaN
else
sixty = mean(zscore(Start:Start+round(60*Fs)));
end

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New 30s Means\', data.info.blockname(1,:)],'thirty', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New 60s Means\', data.info.blockname(1,:)],'sixty', '-ascii');

figure; plot(zscore(Start:(Start+round(30*Fs))));
figure; plot(zscore(Start+(round(30*Fs)+1):Start+round(60*Fs)));
figure; plot(zscore(Start+(round(60*Fs)+1):Start+round(90*Fs)));
figure; plot(zscore(Start+(round(90*Fs)+1):Start+round(120*Fs)));
figure; plot(zscore(Start:Start+round(60*Fs)));

%% Plot and save standard deviation during epocs

stdbefore = std(Before);
stdrotarod = std(Rotarod);
stdafter = std(After);

Stdevs = [stdbefore stdrotarod stdafter];
figure; plot(Stdevs);

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New SDs\', data.info.blockname(1,:)],'Stdevs', '-ascii');

Stdevchange = Stdevs/stdbefore;
stdrotchange = stdrotarod/stdbefore;
figure; plot(Stdevchange);

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\New SD change\', data.info.blockname(1,:)],'Stdevchange', '-ascii');


%% Calculate STDev for 30 and 60s intervals up to 2 minutes on rotarod

if Start+round(30*Fs) > Stop
    thirtySD1 = NaN
else
thirtySD1 = std(normDat(Start:Start+round(30*Fs)));
end

if Start+round(60*Fs) > Stop
    thirtySD2 = NaN
else
thirtySD2 = std(normDat(Start+(round(30*Fs)+1):Start+round(60*Fs)));
end

if Start + round(90*Fs) > Stop
    thirtySD3 = NaN
else
thirtySD3 = std(normDat(Start+(round(60*Fs)+1):Start+round(90*Fs)));
end

if Start + round(120*Fs) > Stop
    thirtySD4 = NaN
else
thirtySD4 = std(normDat(Start+(round(90*Fs)+1):Start+round(120*Fs)));
end

thirtySD = [thirtySD1 thirtySD2 thirtySD3 thirtySD4];
thirtySDchange = thirtySD/stdbefore
figure; plot(thirty)

if Start + round(60*Fs) > Stop
    sixtySD = NaN
else
sixtySD = std(normDat(Start:Start+round(60*Fs)));
sixtySDchange = sixtySD/stdbefore
end

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\30s SDs\', data.info.blockname(1,:)],'thirtySD', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\60s SDs\', data.info.blockname(1,:)],'sixtySD', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\30s SDs change\', data.info.blockname(1,:)],'thirtySDchange', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\2-3 month GCAMP YAC128-FVB\60s SDs change\', data.info.blockname(1,:)],'sixtySDchange', '-ascii');

%% Calculate slope for first 30 sec, 60 sec and entire time on rotarod


