close all
clear variables
clc

%% Import data from TDT files
% Update with current folder path
data = TDTbin2mat(uigetdir);

%% Create time vector from data
% Fs = sampling rate
% tssub = take out first 10 seconds of data

ts = (0:(length(data.streams.x65G.data)-1))/data.streams.x65G.fs;
Fs=data.streams.x65G.fs;
tssub=ts(round(30*Fs):end);

filename = [data.info.blockname]
%% Extract signal and control data and create plot
% Apply lowpass filter 

green=data.streams.x65G.data(round(30*Fs):end);
red=data.streams.x60R.data(round(30*Fs):end);
ctr_green=data.streams.x05G.data(round(30*Fs):end);
ctr_red=data.streams.x05R.data(round(30*Fs):end);
green=double(green);
red=double(red);
ctr_green=double(ctr_green);
ctr_red=double(ctr_red);
ctr_green_filt=lowpassphotometry(ctr_green,Fs,6);
ctr_red_filt=lowpassphotometry(ctr_red,Fs,6); 

plot(tssub,ctr_green_filt)
hold on
plot(tssub,ctr_red_filt)
plot(tssub,green)
hold on
plot(tssub,red)
xlabel('Time(s)','FontSize',14)
ylabel('mV at detector','FontSize',14)
daspect([1 1 1])
hold off
saveas(gcf, 'Raw data figure.jpg')

%% Normalize signal to control channel, plot and save figure
% Add notes to plot

[normGreen] = deltaFF(green,ctr_green_filt);
figure;plot(tssub, normGreen)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('deltaF/F','FontSize',14)
title(filename,'Interpreter','none')
%vline([data.epocs.Note.onset(1,1) data.epocs.Note.onset(2,1) data.epocs.PtAB.onset(1,1) data.epocs.PtAB.onset(2,1)], {'k', 'k', 'g', 'g'});

saveas(gcf, 'Figure GCAMP dFF and Notes.jpg')

hold off

[normRed] = deltaFF(red,ctr_red_filt);
figure;plot(tssub, normRed)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('deltaF/F','FontSize',14)
title(filename,'Interpreter','none')
%vline([data.epocs.Note.onset(1,1) data.epocs.Note.onset(2,1) data.epocs.PtAB.onset(1,1) data.epocs.PtAB.onset(2,1)], {'k', 'k', 'g', 'g'});

saveas(gcf, 'Figure RCAMP dFF and Notes.jpg')

hold off

figure;plot(tssub, normGreen, tssub, normRed);
%vline([data.epocs.Note.onset(1,1) data.epocs.Note.onset(2,1) data.epocs.PtAB.onset(1,1) data.epocs.PtAB.onset(2,1)], {'k', 'k', 'g', 'g'});

%%
Notepoints = (data.epocs.Note.onset.*Fs);
Notepoints = Notepoints - round(10*Fs);
Rotarodpoints = (data.epocs.PtAB.onset.*Fs);
Rotarodpoints = Rotarodpoints - round(10*Fs); 

%Notesoff = (data.epocs.Note.offset.*Fs);
%Notesoff = Notesoff - 10000;
Pickup = round(Notepoints(1,1));
Start = round(Rotarodpoints(1,1));
Stop = round(Rotarodpoints(2,1));
Down = round(Notepoints(2,1));
%End = round(Notesoff(4,1));

%% Calculate New z score using baseline and create plot
Before = normGreen((Pickup-round(300*Fs):(Pickup-round(10*Fs))));
Up = normGreen(Pickup:Start);
Rotarod = normGreen(Start:Stop);
Putdown = normGreen(Stop:Down);
After = normGreen(Down:Down+round(300*Fs));

stdev=std(Before);
zscore=normGreen/stdev;

figure;plot(tssub, zscore)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('Z Score','FontSize',14)
title(filename,'Interpreter','none')
vline([data.epocs.Note.onset(1,1) data.epocs.Note.onset(2,1) data.epocs.PtAB.onset(1,1) data.epocs.PtAB.onset(2,1)], {'k', 'k', 'g', 'g'});

saveas(gcf, 'New Z Score and Notes.jpg')

hold off

%% Downsample data

filt_zscore=lowpassphotometry(zscore,Fs,1);
zscore_down=downsample(filt_zscore,1000)';
ts_down=downsample(tssub,1000)';
%figure; plot(ts_down, zscore_down)

%% Set "start" as 0, plot, save figure and raw data
% Second note should be "start" 
tssub_start = tssub - data.epocs.PtAB.onset(1,1);
note_start = data.epocs.Note.onset - data.epocs.PtAB.onset(1,1);
rotarod_stop = data.epocs.PtAB.onset(2,1) - data.epocs.PtAB.onset(1,1);

figure;plot(tssub_start, zscore)
hold on
xlabel('Time(s)','FontSize',14)
ylabel('Z Score','FontSize',14)
title(filename,'Interpreter','none')
vline([note_start(1,1) note_start(2,1) 0 rotarod_stop(1,1)], {'k', 'k', 'g', 'g'});
daspect([1 0.05 1])
hold off
saveas(gcf, 'Figure Z Score Start Align.jpg');

tssub_start_down=downsample(tssub_start,1000)';

startalign = NaN(300,3);
startalign(1:length(tssub_start_down),1) = tssub_start_down;
startalign(1:length(zscore_down),2) = zscore_down;
startalign(1:length(note_start),3) = note_start;
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New Start Align\', filename],'startalign', '-ascii');

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

Meansnew = [Beforemean Pickupmean Trialmean Putdownmean Aftermean]
figure; plot(Meansnew)
Meanshift = Meansnew - Beforemean
figure; plot(Meanshift)

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New Z Scores\', filename],'Meansnew', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New Norm Z Scores\', filename],'Meanshift', '-ascii');

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

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New 30s Means\', filename],'thirty', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New 60s Means\', filename],'sixty', '-ascii');

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

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New SDs\', filename],'Stdevs', '-ascii');

Stdevchange = Stdevs/stdbefore;
stdrotchange = stdrotarod/stdbefore;
figure; plot(Stdevchange);

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\New SD change\', filename],'Stdevchange', '-ascii');


%% Calculate STDev for 30 and 60s intervals up to 2 minutes on rotarod

if Start+round(30*Fs) > Stop
    thirtySD1 = NaN
else
thirtySD1 = std(normGreen(Start:Start+round(30*Fs)));
end

if Start+round(60*Fs) > Stop
    thirtySD2 = NaN
else
thirtySD2 = std(normGreen(Start+(round(30*Fs)+1):Start+round(60*Fs)));
end

if Start + round(90*Fs) > Stop
    thirtySD3 = NaN
else
thirtySD3 = std(normGreen(Start+(round(60*Fs)+1):Start+round(90*Fs)));
end

if Start + round(120*Fs) > Stop
    thirtySD4 = NaN
else
thirtySD4 = std(normGreen(Start+(round(90*Fs)+1):Start+round(120*Fs)));
end

thirtySD = [thirtySD1 thirtySD2 thirtySD3 thirtySD4];
thirtySDchange = thirtySD/stdbefore
figure; plot(thirty)

if Start + round(60*Fs) > Stop
    sixtySD = NaN
else
sixtySD = std(normGreen(Start:Start+round(60*Fs)));
sixtySDchange = sixtySD/stdbefore
end

save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\30s SDs\', filename],'thirtySD', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\60s SDs\', filename],'sixtySD', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\30s SDs change\', filename],'thirtySDchange', '-ascii');
save(['Z:\Raymond Lab\Ellen\Fiber Photometry\6-7 month GCAMP YAC128\6-7 month rotarod analysis\60s SDs change\', filename],'sixtySDchange', '-ascii');


%% Calculate slope for first 30 sec, 60 sec and entire time on rotarod

%%
openvar('Meansnew');
openvar('stdrotchange');
