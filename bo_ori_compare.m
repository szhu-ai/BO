clear all
close all
clc
%% add the repositories to your path
rootpath='/Users/shushu/Dropbox/npix';
addpath(genpath(fullfile(rootpath,'code')));
testtype='ori';
myresultPath=fullfile(rootpath,'spike',testtype);
myEventDir=fullfile(rootpath,'event_time',testtype); %%% location for event on time 
%% set parameters for the timing counting window for different analysis
sessionidx=5;
flag_exist_celltype=1;
flag_exist_frmat=1;
flag_exist_orifit=1;
flag_exist_MI=1;

p_vr_thresh=0.01;
fr_thresh=3; %% choose firing rate of the optimal condition is at least 1 spikes/sec

%% load basic test information from look up table, please edit this part accordingly
load(fullfile(rootpath,'result','ori_lut.mat'));
recordingDate=ST.recordingDate{sessionidx};
recordingSession=ST.recordingSession{sessionidx};
bankid=ST.bankid(sessionidx);
N_test=ST.Ntest(sessionidx);
eye_tested=ST.eyeID(sessionidx);
mysaveplotPath=fullfile(rootpath,'result',testtype,recordingDate,recordingSession);

depth_maxselection=ST.ymax{sessionidx};
depth_isdeep=ST.isDeep(sessionidx);
depth_relativezero=ST.Zero{sessionidx};
depth_all=[ST.L56{ST.IncludedSession>0};ST.L4c{ST.IncludedSession>0};ST.L4b{ST.IncludedSession>0};ST.L23{ST.IncludedSession>0}];
depth_all_mean_m1=nanmean(depth_all(:,1:3),2);
depth_all_mean_m2=nanmean(depth_all(:,4:5),2);
depth_all(4,[4,5])=nan;
depth_all_mean=nanmean(depth_all,2);
depth_layer_relative=[-depth_all_mean(1);0;depth_all_mean(2);depth_all_mean(2)+depth_all_mean(3);depth_all_mean(2)+depth_all_mean(3)+depth_all_mean(4)];

%% load condition id for all the trials
if  N_test==2
    load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_all1.mat'));   
    trialmat1=load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_trialmat1.mat'));
    trialmat2=load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_trialmat2.mat'));
    trialmat=[trialmat1.trialmat;trialmat2.trialmat];
elseif N_test==1 && bankid==2
    load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_all2.mat'));   
    load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_trialmat2.mat'));
elseif N_test==1 && bankid==1
    load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_all.mat'));   
    load(fullfile(rootpath,'event',recordingDate,recordingSession,'grating_ori_trialmat.mat'));
end
c_condition=4;
c_eye=5;
c_ori=6;
c_sf=7;
N_stim=N_stim*N_test;
N_repetition=N_stim/N_condition;
sf=[0.5,1,2,4];
%% Loading data from kilosort/phy easily  and Computing some useful details about spikes/neurons (like depths)
load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_sp.mat']));
[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
%% load event time
S_event=load(fullfile(myEventDir,recordingDate,recordingSession,['EventOn_',testtype,num2str(bankid),'.mat'])); % a vector of times in seconds of some event to align to 
event_on_time=S_event.event_on_time;
%% compute the depth info for each cluster and sort depth
tempPerClu=findTempForEachClu(sp.clu,sp.spikeTemplates);
tempPerClu(isnan(tempPerClu))=setdiff(1:max(sp.spikeTemplates),unique(sp.spikeTemplates)); 
cluster_templateid=tempPerClu(sp.cids+1); %%% zero-indexed
cluster_depth=templateYpos(cluster_templateid+1);
[cluster.depthsort,cluster.sortidx]=sort(cluster_depth'); 
cluster.depthsorted_id=sp.cids(cluster.sortidx);
cluster.depthsorted_label=sp.cgs(cluster.sortidx);  %% 1 for mu, 2 for su, 3 for noise
cluster.depthsorted_templateid=cluster_templateid(cluster.sortidx);
N_cluster=length(cluster.depthsorted_id);
%%
if flag_exist_celltype==1
    a=load(fullfile(myresultPath,recordingDate,recordingSession,'cell_type_index_shude.mat'));
    b=fieldnames(a);
    cell_type_index=getfield(a,b{1});
    cluster.depthsorted_celltype=cell_type_index(cluster.depthsorted_templateid+1);
else
    cluster.depthsorted_celltype=zeros(1,N_cluster);    
end
N_celltype=numel(unique(cluster.depthsorted_celltype));
if N_celltype==4
%     colorlabel=[55 126 184;77 175 74;154,50,188;228 26 28]/255;
%     colorlabel=[55 126 184;77 175 74;77 175 74;228 26 28]/255;
%     colorlabel=[103,169,207;0 158 115;0 158 115;204 121 167]/255;
%         colorlabel=[103,169,207;86 180 233;86 180 233;230 159 0]/255;
    ccmap=brewermap(11,'PiYG');
    colorlabel_celltype=[67/255,147/255,195/255;ccmap(end-2,:);ccmap(2,:);ccmap(2,:)];
elseif N_celltype==3
    colorlabel_celltype=[0 0 1;0 1 0; 1 0 0];
elseif N_celltype==1
    colorlabel_celltype=[0 0 0];
end
colorlabel_layer=[brewermap(12,'Paired')];
colorlabel_layer=colorlabel_layer([10,4,3,8,2],:);

eyelabel={'L','R','B'};
unitlabel={'MU','SU','N'};
layerline_col=[0.3,0.3,0.3,0.5];
layerline_width=1.5;
mmline_color=[0.7 0.7 0.7 1];
mmline_width=1;
sz_scatter=6;
sz_font=8;

%% compute the firing rate for each trial and for each condition
window_spikecount=[0.1,1]; %% in seconds, relative to stim onset(0)
window_spikecount_bsl=[-0.1,0.1];  %% time window for baseline
psth_window=[-0.1,1.2]; %% time window for plottng psth
psth_binSize=0.01;
binborder_psth = psth_window(1):psth_binSize:psth_window(2);
bincenter_psth = binborder_psth(1:end-1)+psth_binSize/2;
N_Bins = length(binborder_psth)-1; 

if flag_exist_frmat==0
    psth_all=zeros(N_condition,N_Bins,N_cluster);
    psth_std=zeros(N_condition,N_Bins,N_cluster);
    FR_all=zeros(N_stim,N_cluster);
    FR_bsl_all=FR_all;
    for i=1:N_cluster
        idx=0;
        sptime=sp.st(sp.clu==cluster.depthsorted_id(i));
         for id_eye=1:N_eye
            for id_condition=1:N_ori*N_sf
                idx=idx+1;
                event_on_time_currentcond=event_on_time(trialmat(:,c_condition)==id_condition&trialmat(:,c_eye)==id_eye);
                [psth, bincenter_psth, rasterX, rasterY, spikeCounts, cba] = psthAndBA(sptime, event_on_time_currentcond, psth_window, psth_binSize);
                psth_all(idx,:,i) = mean(cba)./psth_binSize;
                psth_std(idx,:,i) = std(cba./psth_binSize,1);
            end
        end       
        S1=arrayfun(@(x) length(sptime(sptime<=x+window_spikecount(2)&sptime>=x+window_spikecount(1))), event_on_time,'UniformOutput',false);
        FR_all(:,i)=cell2mat(S1)./(window_spikecount(2)-window_spikecount(1));
        S0=arrayfun(@(x) length(sptime(sptime<=x+window_spikecount_bsl(2)&sptime>=x+window_spikecount_bsl(1))), event_on_time,'UniformOutput',false);
        FR_bsl_all(:,i)=cell2mat(S0)./(window_spikecount_bsl(2)-window_spikecount_bsl(1));
    end 
    [~,sortidx]=sort(trialmat(:,c_condition)+(trialmat(:,c_eye)-1)*N_ori*N_sf);
    FR_all_sorted=FR_all(sortidx,:);
    FR_bsl_all_sorted=FR_bsl_all(sortidx,:);
    aa=reshape(FR_all_sorted,N_repetition,N_condition,N_cluster);
    FR_mat=squeeze(mean(aa,1));
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_psth_all.mat']),'psth_all');   
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_psth_std.mat']),'psth_std');        
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRall.mat']),'FR_all');
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRall_sorted.mat']),'FR_all_sorted');
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRbsl_all.mat']),'FR_bsl_all');    
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRall_bsl_sorted.mat']),'FR_bsl_all_sorted');   
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRcond.mat']),'FR_mat');
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_cluster.mat']),'cluster');
else
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_psth_all.mat']));   
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_psth_std.mat']));     
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRall.mat']));
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRall_sorted.mat']));
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRbsl_all.mat']));
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRall_bsl_sorted.mat'])); 
    load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_FRcond.mat']));    
end
%% determin the relative depth
cell_layer_idx=nan(length(cluster.depthsort),1);
depth_scaled_temp=cluster.depthsort;
if depth_isdeep==0
    depth_L56=[depth_relativezero-ST.L56{sessionidx},depth_relativezero];
    depth_L4c=[depth_relativezero,depth_relativezero+ST.L4c{sessionidx}];
    depth_L4b=[depth_L4c(2),depth_L4c(2)+ST.L4b{sessionidx}];
    depth_L23=[depth_L4b(2),depth_L4b(2)+ST.L23{sessionidx}];
    depth_WM=[0,depth_L56(1)];
elseif depth_isdeep==1
    depth_L56=[depth_relativezero,depth_relativezero+ST.L56{sessionidx}];
    depth_L4c=[depth_relativezero-ST.L4c{sessionidx},depth_relativezero];
    depth_L4b=[depth_L4c(1)-ST.L4b{sessionidx},depth_L4c(1)];
    depth_L23=[depth_L4b(1)-ST.L23{sessionidx},depth_L4b(1)];
%     depth_L4b(1)=max(depth_L4b(1),0);
%     if depth_L4b(1)>0
%         depth_L23(1)=max(depth_L23(1),0);
%     end
    depth_WM=[depth_L56(2),depth_maxselection];
end
depth_mat=[depth_L56;depth_L4c;depth_L4b;depth_L23;depth_WM];

N_layer=size(depth_mat,1);
for id_layer=1:N_layer
    cell_layer_idx(cluster.depthsort>=depth_mat(id_layer,1)&cluster.depthsort<depth_mat(id_layer,2))=id_layer;               
end

%% paired t test for visual responseviness
[FR_ttest_h,FR_ttest_p]=ttest(FR_all_sorted,FR_bsl_all_sorted,'Alpha',p_vr_thresh,'Tail','right');
Includedidx=FR_ttest_h>0&max(FR_dir_prefsf)>fr_thresh&cluster.depthsorted_celltype>0;