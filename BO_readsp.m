
%% add the repositories to your path
clear all
clc

%% set paths for where to find your data
sessionidx=4;
rootpath='G:\npix\';
testtype='BO';
myresultPath=fullfile(rootpath,'spike',testtype);
addpath(genpath(fullfile(rootpath,'code')));
load(fullfile(rootpath,'code','bo_lut.mat'));
recordingDate=ST.recordingDate{sessionidx};
recordingSession=ST.recordingSession{sessionidx};
bankid=ST.bankid(sessionidx);
id_eyelable=ST.eyeID(sessionidx);

myKsDir = 'G:\Quicky_kilo2_bo\';
myEventDir=fullfile(rootpath,'event_time',testtype); %%% location for event on time 

%% Loading data from kilosort/phy easily
sp = loadKSdir(fullfile(myKsDir,recordingDate,recordingSession)); % % sp.st are spike times in seconds; sp.clu are cluster identities, spikes from clusters labeled "noise" have already been omitted
% sp = loadKSdir(fullfile(myKsDir)); % % sp.st are spike times in seconds; sp.clu are cluster identities, spikes from clusters labeled "noise" have already been omitted
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_sp.mat']),'sp');

%% load synchronization data
flag_read_event=0;
if flag_read_event==1
    syncChanIndex = sp.n_channels_dat;
    nChansInFile = sp.n_channels_dat;  % neuropixels phase3a, from spikeGLX
%     syncDat = extractSyncChannel([myKsDir,recordingDate,'\',recordingSession], nChansInFile, syncChanIndex);
    syncDat = extractSyncChannel(['G:\npix\raw\2018_12_06\Data\ori1'], nChansInFile, syncChanIndex);

    de2bi=1;
    apFs=sp.sample_rate;
%     apFs=2500;
    myEventTimes = spikeGLXdigitalParse(syncDat, apFs);
    % - eventTimes{1} contains the sync events from digital channel 1, as three cells: 
    % - eventTimes{1}{1} is the times of all events
    % - eventTimes{1}{2} is the times the digital bit went from off to on
    % - eventTimes{1}{2} is the times the digital bit went from on to off
    event_on=cell(1,2);
    for DIchan=[1 2]
        event_on_all=myEventTimes{1,DIchan}{1,2};
        for eventidx=2:length(event_on_all)
            if event_on_all(eventidx)-event_on_all(eventidx-1)<=0.6     % manully select a threshold for trigger interval and determine the total number to match experimental one.
                event_on_all(eventidx)=event_on_all(eventidx-1);
            end
        end
        devent_on=[event_on_all(1);diff(event_on_all)];
        event_on{1,DIchan}=event_on_all(devent_on>0);
    end
    event_on_time=event_on{1,1}(1:end);
    save(fullfile(myEventDir,recordingDate,recordingSession,['EventOn_',testtype,num2str(bankid),'.mat']),'event_on_time'); % a vector of times in seconds of some event to align to 
    save(fullfile(myEventDir,recordingDate,recordingSession,['EventOn_',testtype,num2str(bankid),'.mat']),'myEventTimes'); % a vector of times in seconds of some event to align to 
end

%%
for i=1:10
    aa=max(abs(sp.temps(i,:,:)),[],2);
    [~,ab]=max(aa);
    figure;
    plot(sp.temps(i,:,ab))
end