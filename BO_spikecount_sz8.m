clear all
close all
clc
for sessionidx=13% [16 15 14 13 11] %[5 6 15 16 10 13 11 7 14]
    clearvars -except sessionidx
    clc
    fprintf('Sessionidx %d.\n',sessionidx)
    %% add the repositories to your path
    rootpath='/Users/shushu/Dropbox/npix';
    addpath(genpath(fullfile(rootpath,'code')));
    addpath(genpath(fullfile(rootpath,'code_new')));
    testtype='BO';
    myresultPath=fullfile(rootpath,'spike',testtype);
    myEventDir=fullfile(rootpath,'event_time',testtype); %%% location for event on time 
    %% set parameters for the timing counting window for different analysis
    flag_exist_frmat=1;
    flag_exist_celltype=0;
    flag_plot_ori_corr=0;
    flag_plot_raster=0;
    flag_exis_layer=1;
    flag_exist_decoding=1;
    flag_plot_decoding=0;
    flag_exist_erp=1;
    flag_plot_csd=1;
    flag_exist_corr=1;
    flag_plot_ccg=0;

    p_vr_thresh=0.05;
    fr_thresh=0.5; %% choose firing rate of the optimal condition is at least 1 spikes/sec

    window_psth=[-0.3,1]; %% time window for plottng psth
    psth_binSize=0.001;
    filter_window=0.05;
    binborder_psth = window_psth(1):psth_binSize:window_psth(2);
    bincenter_psth = binborder_psth(1:end-1)+psth_binSize/2;
    N_Bins = length(binborder_psth)-1; 
    
    window_bsl_psth=[-0.3,0]; %% time window for plottng psth
    binborder_bsl_psth = window_bsl_psth(1):psth_binSize:window_bsl_psth(2);
    N_Bins_bsl = length(binborder_bsl_psth)-1; 
    
    window_spikecount=[0.05,0.5]; %% in seconds, relative to stim onset(0)
    window_spikecount_idx=find(binborder_psth>=window_spikecount(1),1,'first') :1: find(binborder_psth<=window_spikecount(2),1,'last')-1;

    colorlabel_layer=brewermap(12,'Paired');
    colorlabel_layer=colorlabel_layer([10,4,3,8,2],:);
    %% load basic test information from look up table, please edit this part accordingly
%     load(fullfile(rootpath,'result','bo_lut.mat'));
    load(fullfile(rootpath,'result','bo_lut_based_on_ori_depth.mat'));
    ori_lut=load(fullfile(rootpath,'result','ori_lut_based_on_bo_depth.mat'));

    recordingDate=ST.recordingDate{sessionidx};
    recordingSession=ST.recordingSession{sessionidx};
    bankid=ST.bankid(sessionidx);
    N_test=ST.Ntest(sessionidx);
    eye_tested=ST.eyeID(sessionidx);
    mysaveplotPath=fullfile(rootpath,'result',testtype,recordingDate,recordingSession);
    depth_maxselection=ST.ymax{sessionidx};
    depth_isdeep=ST.isDeep(sessionidx);
    %% calculate depth based on orientation session
    depth_relativezero=ori_lut.ST.Zero{sessionidx};
    if depth_isdeep==0
        depth_L56=[depth_relativezero-ori_lut.ST.L56{sessionidx},depth_relativezero];
        depth_L4c=[depth_relativezero,depth_relativezero+ori_lut.ST.L4c{sessionidx}];
        depth_L4b=[depth_L4c(2),depth_L4c(2)+ori_lut.ST.L4b{sessionidx}];
        depth_L23=[depth_L4b(2),depth_L4b(2)+ori_lut.ST.L23{sessionidx}];
        depth_WM=[0,depth_L56(1)];

        depth_full_L6=[depth_relativezero-ori_lut.ST.L56{sessionidx},depth_relativezero-ori_lut.ST.L56{sessionidx}/2];
        depth_full_L5=[depth_full_L6(2),depth_relativezero];
        depth_full_L4cb=[depth_full_L5(2),depth_full_L5(2)+ori_lut.ST.L4c{sessionidx}/2];
        depth_full_L4ca=[depth_full_L4cb(2),depth_full_L4cb(2)+ori_lut.ST.L4c{sessionidx}/2];
        depth_full_L4b=[depth_full_L4ca(2),depth_full_L4ca(2)+ori_lut.ST.L4b{sessionidx}];
        depth_full_L3=[depth_full_L4b(2),depth_full_L4b(2)+ori_lut.ST.L23{sessionidx}/2];
        depth_full_L2=[depth_full_L3(2),depth_full_L3(2)+ori_lut.ST.L23{sessionidx}/2];
    elseif depth_isdeep==1
        depth_L56=[depth_relativezero,depth_relativezero+ori_lut.ST.L56{sessionidx}];
        depth_L4c=[depth_relativezero-ori_lut.ST.L4c{sessionidx},depth_relativezero];
        depth_L4b=[depth_L4c(1)-ori_lut.ST.L4b{sessionidx},depth_L4c(1)];
        depth_L23=[depth_L4b(1)-ori_lut.ST.L23{sessionidx},depth_L4b(1)];
        depth_L4b(1)=max(depth_L4b(1),0);
        if depth_L4b(1)>0
            depth_L23(1)=max(depth_L23(1),0);
        end
        depth_WM=[depth_L56(2),depth_maxselection];

        depth_full_L6=[depth_relativezero+ori_lut.ST.L56{sessionidx}/2,depth_relativezero+ori_lut.ST.L56{sessionidx}];
        depth_full_L5=[depth_relativezero,depth_full_L6(1)];
        depth_full_L4cb=[depth_full_L5(1)-ori_lut.ST.L4c{sessionidx}/2,depth_full_L5(1)];
        depth_full_L4ca=[depth_full_L4cb(1)-ori_lut.ST.L4c{sessionidx}/2,depth_full_L4cb(1)];
        depth_full_L4b=[depth_full_L4ca(1)-ori_lut.ST.L4b{sessionidx},depth_full_L4ca(1)];
        depth_full_L3=[depth_full_L4b(1)-ori_lut.ST.L23{sessionidx},depth_full_L4b(1)];
        depth_full_L4b(1)=max(depth_full_L4b(1),0);
        if depth_full_L4b(1)>0
            depth_full_L3(1)=max(depth_full_L3(1),0);
        end
        depth_full_L2=[-2,-1];
    end
    depth_mat_ori=[depth_L56;depth_L4c;depth_L4b;depth_L23;depth_WM];
    depth_full_mat_ori=[depth_full_L6;depth_full_L5;depth_full_L4cb;depth_full_L4ca;depth_full_L4b;depth_full_L3;depth_full_L2;depth_WM];
    %% calculate depth based on bo session
    depth_relativezero=ST.Zero{sessionidx};
    if depth_isdeep==0
        depth_L56=[depth_relativezero-ST.L56{sessionidx},depth_relativezero];
        depth_L4c=[depth_relativezero,depth_relativezero+ST.L4c{sessionidx}];
        depth_L4b=[depth_L4c(2),depth_L4c(2)+ST.L4b{sessionidx}];
        depth_L23=[depth_L4b(2),depth_L4b(2)+ST.L23{sessionidx}];
        depth_WM=[0,depth_L56(1)];

        depth_full_L6=[depth_relativezero-ST.L56{sessionidx},depth_relativezero-ST.L56{sessionidx}/2];
        depth_full_L5=[depth_full_L6(2),depth_relativezero];
        depth_full_L4cb=[depth_full_L5(2),depth_full_L5(2)+ST.L4c{sessionidx}/2];
        depth_full_L4ca=[depth_full_L4cb(2),depth_full_L4cb(2)+ST.L4c{sessionidx}/2];
        depth_full_L4b=[depth_full_L4ca(2),depth_full_L4ca(2)+ST.L4b{sessionidx}];
        depth_full_L3=[depth_full_L4b(2),depth_full_L4b(2)+ST.L23{sessionidx}/2];
        depth_full_L2=[depth_full_L3(2),depth_full_L3(2)+ST.L23{sessionidx}/2];
    elseif depth_isdeep==1
        depth_L56=[depth_relativezero,depth_relativezero+ST.L56{sessionidx}];
        depth_L4c=[depth_relativezero-ST.L4c{sessionidx},depth_relativezero];
        depth_L4b=[depth_L4c(1)-ST.L4b{sessionidx},depth_L4c(1)];
        depth_L23=[depth_L4b(1)-ST.L23{sessionidx},depth_L4b(1)];
        depth_L4b(1)=max(depth_L4b(1),0);
        if depth_L4b(1)>0
            depth_L23(1)=max(depth_L23(1),0);
        end
        depth_WM=[depth_L56(2),depth_maxselection];

        depth_full_L6=[depth_relativezero+ST.L56{sessionidx}/2,depth_relativezero+ST.L56{sessionidx}];
        depth_full_L5=[depth_relativezero,depth_full_L6(1)];
        depth_full_L4cb=[depth_full_L5(1)-ST.L4c{sessionidx}/2,depth_full_L5(1)];
        depth_full_L4ca=[depth_full_L4cb(1)-ST.L4c{sessionidx}/2,depth_full_L4cb(1)];
        depth_full_L4b=[depth_full_L4ca(1)-ST.L4b{sessionidx},depth_full_L4ca(1)];
        depth_full_L3=[depth_full_L4b(1)-ST.L23{sessionidx},depth_full_L4b(1)];
        depth_full_L4b(1)=max(depth_full_L4b(1),0);
        if depth_full_L4b(1)>0
            depth_full_L3(1)=max(depth_full_L3(1),0);
        end
        depth_full_L2=[-2,-1];
    end
    depth_mat=[depth_L56;depth_L4c;depth_L4b;depth_L23;depth_WM];
    depth_full_mat=[depth_full_L6;depth_full_L5;depth_full_L4cb;depth_full_L4ca;depth_full_L4b;depth_full_L3;depth_full_L2;depth_WM];
    %% load condition id for all the trials
    if  N_test==2
        load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_all1.mat'));   
        trialmat1=load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_trialmat1.mat'));
        trialmat2=load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_trialmat2.mat'));
        trialmat=[trialmat1.trialmat;trialmat2.trialmat];
    elseif N_test==1 && bankid==2
        load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_all2.mat'));   
        load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_trialmat2.mat'));
    elseif N_test==1 && bankid==1
        load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_all.mat'));   
        load(fullfile(rootpath,'event',recordingDate,recordingSession,'bo_trialmat.mat'));
    end
    % c_eye=5;
    c_ori=6;
    c_side=7;
    c_lc=8;
    c_sz=9;
    trialmat(trialmat(:,c_side)==2,c_lc)=3-trialmat(trialmat(:,c_side)==2,c_lc);
    id_condition=trialmat(:,c_sz)+(trialmat(:,c_side)-1)*N_sz+(trialmat(:,c_lc)-1)*N_sz*N_lc+(trialmat(:,c_ori)-1)*N_sz*N_lc*N_side;
    N_stim=N_stim*N_test;
    N_repetition=size(trialmat,1)*N_test/N_condition;
    N_ori=4;
    %% Loading data from kilosort/phy easily  and Computing some useful details about spikes/neurons (like depths)
%     sp = loadKSdir(fullfile(rootpath,'raw','BO','kilosort3_old',recordingDate,recordingSession,'kilosort3_curated'));
    sp = loadKSdir(fullfile('/Users/shushu/kilosort/bo',recordingDate,recordingSession,'kilosort3'));
    save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_sp.mat']),'sp','-v7.3');   

%     [spikeAmps, spikeDepths, templateXpos,templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
%     templatePositionsAmplitudes_shude(sp.temps, sp.winv, sp.xcoords,sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
    [~, ~, templateXpos,templateYpos, tempAmps, ~, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes_shude(sp.temps, sp.winv, sp.xcoords,sp.ycoords, sp.spikeTemplates, sp.tempScalingAmps);
    if size(sp.cids,1)>size(sp.cids,2)
        sp.cids=sp.cids';
    end
    if size(sp.cgs,1)>size(sp.cgs,2)
        sp.cgs=sp.cgs';
    end
    %% load event time
    S_event=load(fullfile(myEventDir,recordingDate,recordingSession,['EventOn_',testtype,num2str(bankid),'.mat'])); % a vector of times in seconds of some event to align to 
    event_on_time=S_event.event_on_time;
    if sessionidx==7
        event_on_time=event_on_time(1:length(id_condition));
    end
    %% compute the depth info for each cluster and sort in accending manner, depth 0 is the tip of the probe
    tempPerClu=findTempForEachClu(sp.clu,sp.spikeTemplates);
%     tempPerClu(isnan(tempPerClu))=setdiff(1:max(sp.spikeTemplates),unique(sp.spikeTemplates)); 
    cluster_templateid=tempPerClu(sp.cids+1); %%% zero-indexed
    cluster_xpos=templateXpos(cluster_templateid+1);   % cluster_depth=templateYpos;
    cluster_xpos=cluster_xpos';
    cluster_depth=templateYpos(cluster_templateid+1);   % cluster_depth=templateYpos;
    [cluster.depthsort,cluster.sortidx]=sort(cluster_depth'); 
    cluster.depthsorted_xpos=cluster_xpos(cluster.sortidx);
    cluster.depthsorted_id=sp.cids(cluster.sortidx);
    cluster.depthsorted_label=sp.cgs(cluster.sortidx);  %% 1 for su, 2 for mu, 3 for noise
    cluster.depthsorted_templateid=cluster_templateid(cluster.sortidx);
    cluster_amps=tempAmps(cluster_templateid+1);
    cluster.depthsorted_amps=cluster_amps(cluster.sortidx);
    N_cluster=length(cluster.depthsorted_id);
    cluster_peakWF=tempPeakWF(cluster_templateid+1,:);
    cluster.depthsorted_peakWF=cluster_peakWF(cluster.sortidx,:);
    %%
    if flag_exist_celltype==1
%         load(fullfile(myresultPath,recordingDate,recordingSession,'cell_type_index_shude.mat'));
%         b=fieldnames(a);
%         cell_type_index=getfield(a,b{1});
%         cluster.depthsorted_celltype=cell_type_index(cluster.depthsorted_templateid+1);
    else
        BO_waveform_V1;
    end
    cluster.depthsorted_celltype=celltype(cluster.depthsorted_templateid+1);
    N_celltype=numel(unique(cluster.depthsorted_celltype));
    if N_celltype==4
        waveform_colorlabel=[55 126 184;77 175 74;154,50,188;228 26 28]/255;
    elseif N_celltype==3
        waveform_colorlabel=[0 0 1;0 1 0; 1 0 0];
    elseif N_celltype==1
        waveform_colorlabel=[0 0 0];
    end
%%
    N_layer=size(depth_mat,1);
    cell_layer_idx=ones(1,length(cluster.depthsort))*-99;
    cell_layer_idx_ori=ones(1,length(cluster.depthsort))*-99;
    for id_layer=1:N_layer
        cell_layer_idx(cluster.depthsort>=depth_mat(id_layer,1)&cluster.depthsort<depth_mat(id_layer,2))=id_layer;  
        cell_layer_idx_ori(cluster.depthsort>=depth_mat_ori(id_layer,1)&cluster.depthsort<depth_mat_ori(id_layer,2))=id_layer;          
    end

    N_layer_full=size(depth_full_mat,1);
    cell_layer_full_idx=ones(1,length(cluster.depthsort))*-99;
    cell_layer_full_idx_ori=ones(1,length(cluster.depthsort))*-99;
    for id_layer=1:N_layer_full
        cell_layer_full_idx(cluster.depthsort>=depth_full_mat(id_layer,1)&cluster.depthsort<depth_full_mat(id_layer,2))=id_layer;          
        cell_layer_full_idx_ori(cluster.depthsort>=depth_full_mat_ori(id_layer,1)&cluster.depthsort<depth_full_mat_ori(id_layer,2))=id_layer;          
    end
    %% compute the firing rate for each trial, saved in FR_all, and for each condition saved in FR_currentcond
    smoothSize = 10; % in msec, stdev of gaussian smoothing filter
    gw = gausswin(round(smoothSize*6),3);
    smWin = gw./sum(gw);    
    if flag_exist_frmat==0
        data=nan(length(event_on_time),N_Bins,N_cluster); 
        data_bsl=nan(length(event_on_time),N_Bins_bsl,N_cluster);         
        sptime_all=cell(N_cluster,1);
        N_spike_raw=zeros(N_cluster,1);
        for id_cluster=1:N_cluster
            sptime=sp.st(sp.clu==cluster.depthsorted_id(id_cluster));
            sptime((diff(sptime)<=0.00016))=[];
            sptime_all{id_cluster}=sptime;
            N_spike_raw(id_cluster)=length(sptime_all{id_cluster});
        end
        for id_cluster=1:N_cluster
            sptime=sptime_all{id_cluster};
            cnt=find(cluster.depthsort-cluster.depthsort(id_cluster)<=50 & cluster.depthsort-cluster.depthsort(id_cluster)>=0);
            if ~isempty(cnt)
                for j=1:length(cnt)
                    clu2_id=cnt(j);
                    if cluster.depthsorted_id(clu2_id)>cluster.depthsorted_id(id_cluster)
                       sptime2=sort(sptime_all{clu2_id});
                       sptime3=sort([sort(sptime);sptime2]);
                       doublecounted_idx1=find(diff(sptime3)<=0.00016);
                       doublecounted=sptime3([doublecounted_idx1,doublecounted_idx1+1]);
                       if cluster.depthsorted_amps(id_cluster)<cluster.depthsorted_amps(clu2_id)
                          sptime(ismember(sptime,doublecounted))=[];
                          sptime_all{id_cluster}=sptime;
                       else
                          sptime2(ismember(sptime2,doublecounted))=[];
                          sptime_all{clu2_id}=sptime2;
                       end
                    end
                end
            end
        end
        for id_cluster=1:N_cluster
            id_cluster
            sptime=sptime_all{id_cluster};
            [psth, bincenter_psth, rasterX, rasterY, spikeCounts, data(:,:,id_cluster)] = psthAndBA(sptime, event_on_time, window_psth, psth_binSize); 
            [~, ~, ~, ~, ~, data_bsl(:,:,id_cluster)] = psthAndBA(sptime, event_on_time, window_bsl_psth, psth_binSize);       
        end 
        save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_data.mat']),'data','-v7.3');   
        save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_data_bsl.mat']),'data_bsl','-v7.3');   
%         save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_psth_cond_all.mat']),'psth_cond_all','-v7.3');           
    else
        load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_data.mat']),'data'); 
        load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_data_bsl.mat']),'data_bsl');       
%         load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_psth_cond_all.mat']),'psth_cond_all');          
    end
    %% determine the visual responseviness and exlcude neurons
    %%%% for each orientation condition. At least 1 out of 4 conditons (2 lc by 2 side)
    %%%% 1)  has firing rate at least 1spikes/sec
    %%%% 2)  evoked FR is different from BG in all conditions (ranksum, bonferroni corrected)
    FR_cond=nan(N_condition,N_cluster);
    FR_mat=nan(N_condition,N_repetition,N_cluster);
    FR_mat_bsl=nan(N_condition,N_repetition,N_cluster);
    N_outlier=nan(N_condition,N_cluster);
    psth_cond_all=nan(N_condition,N_Bins,N_cluster);
    psth_cond_all_zscore=nan(N_condition,N_Bins,N_cluster);
    for id_cond=1:N_condition
%         id_cond_selected=find(id_condition==id_cond);
%         data_i=mean(data(id_cond_selected,window_spikecount_idx,:),2);
%         data_i_zscore=zscore(data_i);
%         id_outlier=data_i_zscore>=4;
        for id_cluster=1:N_cluster
%             if ~isempty(find(id_outlier(:,:,id_cluster)==1))
%                 N_outlier(id_cond,id_cluster)=length(find(id_outlier(:,:,id_cluster)==1));
%                 data(id_cond_selected(id_outlier(:,:,id_cluster)==1),:,id_cluster)=mean(data(id_cond_selected(id_outlier(:,:,id_cluster)==0),:,id_cluster),1);
%                 data_bsl(id_cond_selected(id_outlier(:,:,id_cluster)==1),:,id_cluster)=mean(data_bsl(id_cond_selected(id_outlier(:,:,id_cluster)==0),:,id_cluster),1);                                                     
%             end
            temp=squeeze(mean(data(id_condition==id_cond,:,id_cluster),1,'omitnan'));
            psth_cond_all(id_cond,:,id_cluster)=conv2(smWin,1,temp', 'same')'./psth_binSize; 
        end
        FR_mat(id_cond,:,:)=squeeze(mean(data(id_condition==id_cond,window_spikecount_idx,:),2,'omitnan'))./psth_binSize;
        FR_mat_bsl(id_cond,:,:)=squeeze(mean(data_bsl(id_condition==id_cond,:,:),2,'omitnan'))./psth_binSize;
        temp=squeeze(mean(data(id_condition==id_cond,window_spikecount_idx,:),[1,2],'omitnan'))./psth_binSize;
        FR_cond(id_cond,:)=temp';              
    end
    FR_bsl_mat=squeeze(mean(data_bsl,2,'omitnan'))./psth_binSize;
    FR_mat_reshape=reshape(FR_mat,8,4,N_repetition,N_cluster);
    FR_mat_reshape_sz8=FR_mat_reshape([2,4,6,8],:,:,:);
    FR_mat_reshape_sz4=FR_mat_reshape([1,3,5,7],:,:,:);
    FR_mat_bsl_reshape=reshape(FR_mat_bsl,8,4,N_repetition,N_cluster);
    FR_mat_bsl_reshape_sz8=FR_mat_bsl_reshape([2,4,6,8],:,:,:);
    [FR_ttest_h,FR_ttest_p]=ttest(squeeze(nanmean(data(mod(id_condition,2)==0,window_spikecount_idx,:),2)),squeeze(nanmean(data_bsl(mod(id_condition,2)==0,:,:),2)),'Tail','right');
%     psth_max=squeeze(max(psth_cond_all(:,1:300,:),[],[1 2]));
    Includedidx_all=FR_ttest_p<p_vr_thresh&cell_layer_idx>0 & max(FR_cond(2:2:32,:))>fr_thresh; % | psth_max'>15);
%     Includedidx_all=FR_ttest_p<0.01&FR_ttest_p<p_vr_thresh&cell_layer_idx>0;% & max(FR_cond(2:2:32,:))<=fr_thresh; % | psth_max'>15);

    Includedidx_all=Includedidx_all';
%     vis_p=nan(4,N_cluster,N_ori);
%     vis_fr=nan(4,N_cluster,N_ori);
%     for id_ori=1:4
%         FR_mat_reshape_sz8_i=squeeze(FR_mat_reshape_sz8(:,id_ori,:,:));
%         FR_mat_bsl_reshape_sz8_i=squeeze(FR_mat_bsl_reshape_sz8(:,id_ori,:,:));
%         for idx=1:4
%             for id_cluster=1:N_cluster
%                 vis_p(idx,id_cluster,id_ori)=ranksum(FR_mat_reshape_sz8_i(idx,:,id_cluster),FR_mat_bsl_reshape_sz8_i(idx,:,id_cluster),'Tail','right');
%                 vis_fr(idx,id_cluster,id_ori)=nanmean(FR_mat_reshape_sz8_i(idx,:,id_cluster));
%             end
%         end
%     end
%     Includedidx=squeeze(min(vis_p))<0.05/16&squeeze(max(vis_fr))>0&cell_layer_idx'>0;
%     Includedidx_all=any(Includedidx,2);

    Ncell_included_all=sum(Includedidx_all)
    FRmax=max(FR_cond(2:2:32,Includedidx_all));FRmax=FRmax';
%%
    Ncell_layer=zeros(N_layer,1);
    for id_layer=1:N_layer
        Ncell_layer(id_layer)=length(find(cell_layer_idx==id_layer&Includedidx_all'));             
    end
 %%

%     data2=data;
%     clear data
%     data.spiketrain=data2(:,301:end,:);
%     data.psth_window=[0,1];
%     data.cluster_id=cluster.depthsorted_id';
%     data.cluster_ypos=cluster.depthsort;   %% optional
%     data.layerID=cell_layer_idx;
%     data.Includedidx=Includedidx_all';
%     save(['/Users/shushu/Dropbox/npix/code_new/CCG_V1/data/data_BO_',recordingDate,recordingSession,'.mat'],'data','-v7.3');
%     % data.Includedidx=Includedidx;
%%
    cluster.depthsorted_id_includedall=cluster.depthsorted_id(:,Includedidx_all);
    cluster.depthsorted_label_includedall=cluster.depthsorted_label(:,Includedidx_all);
    cluster.depthsort_includedall=cluster.depthsort(1,Includedidx_all);
    cluster.depthsorted_celltype_includedall=cluster.depthsorted_celltype(1,Includedidx_all);
    cluster.depthsorted_xpos_includedall=cluster.depthsorted_xpos(1,Includedidx_all);
    cluster.depthsorted_peakWF_includedall=cluster.depthsorted_peakWF(Includedidx_all,:);
%     cell_layer_idx_border_included_i=cell(1,4);
%     for id_ori=1:4
%         cell_layer_idx_included_i=cell_layer_idx(Includedidx(:,id_ori));
%         cell_layer_idx_border_included_i{id_ori}=find(diff(cell_layer_idx_included_i));
%     end
    cluster.depthsorted_celllayer_includedall_ori=cell_layer_idx_ori(Includedidx_all);
    cluster.depthsorted_celllayer_full_includedall_ori=cell_layer_full_idx_ori(Includedidx_all);

    cluster.depthsorted_celllayer_includedall=cell_layer_idx(Includedidx_all);
    cell_layer_idx_includedall=cell_layer_idx(Includedidx_all);  
    cell_layer_idx_border_included=find(diff(cluster.depthsorted_celllayer_includedall));

    cluster.depthsorted_celllayer_full_includedall=cell_layer_full_idx(Includedidx_all);
    cell_layer_full_idx_includedall=cell_layer_full_idx(Includedidx_all);   
    cell_layer_full_idx_border_included=find(diff(cluster.depthsorted_celllayer_full_includedall));
    %%
    if flag_exist_decoding==0
        BO_decoding_V2
    end
    %%
    if flag_exist_corr==0
        BO_ccg_v1;
    end
end
%% plot waveform
figure('Color',[1 1 1],'WindowState', 'normal');
for id_celltype=[0,3,2,1]
    plot(repmat(cluster.depthsorted_xpos_includedall(cluster.depthsorted_celltype_includedall==id_celltype),N_waveform_timepoints,1)+dx*repmat(time_wf',1,length(find(cluster.depthsorted_celltype_includedall==id_celltype))),repmat(cluster.depthsort_includedall(cluster.depthsorted_celltype_includedall==id_celltype),N_waveform_timepoints,1)+dy*cluster.depthsorted_peakWF_includedall(cluster.depthsorted_celltype_includedall==id_celltype',:)','Color',waveform_colorlabel(id_celltype+1,:),'LineWidth',2.5);
    hold on
end
layerline_width=2;
for id_layer=1:5
    plot([10,70],[depth_layer_upper(id_layer),depth_layer_upper(id_layer)],'Color',[0 0 0 0.5],'LineWidth',layerline_width,'LineStyle','--')
    hold on
end
ax=gca;
ax.XLabel.String='X-location (\mum) ';
ax.YLabel.String='Y-location (\mum) ';
ax.YDir='Normal';
ax.TickDir = 'out';
ax.PlotBoxAspectRatio = [1 1.5 1];
ax.Box='off';
% ax.YLim=[0,depth_maxselection];
% print('-vector','-dpdf', [mysaveplotPath,'\',recordingDate,recordingSession,'_waveform_',testtype,'.pdf'], '-r0');  
%%
if flag_plot_raster==1
    BO_raster_plot_v1
end
if flag_exist_erp==1
    load(fullfile(myresultPath,recordingDate,recordingSession,'ERP_bo1.mat'));   
end
if flag_plot_csd==1
    csd_plot
end
%%
if flag_plot_decoding==1
    BO_decoding_plot_v1;
end
%% run anova with factor "side" and "LC", only use data from size 8
p_sig=0.01;
p_anova=nan(N_cluster,N_ori,3);
p_anova_side=nan(N_cluster,N_ori,2);

BO_index_FR=nan(N_cluster,N_ori,4);

BO_index=nan(N_cluster,N_ori,3);
LC_index=nan(N_cluster,N_ori,3);
BO_index_normbymax=nan(N_cluster,N_ori,3);
LC_index_normbymax=nan(N_cluster,N_ori,3);

FR_cond_sz8=FR_cond(2:2:32,:);
for id_ori=1:N_ori
    id_trial_selected=trialmat(:,c_ori)==id_ori&trialmat(:,c_sz)==2;
    trialmat_i=trialmat(id_trial_selected,:); 
    for id_cluster=1:N_cluster
        if Includedidx_all(id_cluster)
            FR_i=squeeze(nanmean(data(id_trial_selected,window_spikecount_idx,id_cluster),2))./psth_binSize;
            [p_anova(id_cluster,id_ori,:),tbl]=anovan(FR_i,{trialmat_i(:,c_side), trialmat_i(:,c_lc)},'model','interaction','varnames',{'side','lc'},'display','off');
            [p_anova_side(id_cluster,id_ori,1),tbl]=anova1(FR_i(trialmat_i(:,c_lc)==1),trialmat_i(trialmat_i(:,c_lc)==1,c_side),'off');
            [p_anova_side(id_cluster,id_ori,2),tbl]=anova1(FR_i(trialmat_i(:,c_lc)==2),trialmat_i(trialmat_i(:,c_lc)==2,c_side),'off');
            FR_i_side1_lc1=nanmean(FR_i(trialmat_i(:,c_side)==1&trialmat_i(:,c_lc)==1));
            FR_i_side2_lc1=nanmean(FR_i(trialmat_i(:,c_side)==2&trialmat_i(:,c_lc)==1));
            FR_i_side1_lc2=nanmean(FR_i(trialmat_i(:,c_side)==1&trialmat_i(:,c_lc)==2));
            FR_i_side2_lc2=nanmean(FR_i(trialmat_i(:,c_side)==2&trialmat_i(:,c_lc)==2));
            temp_max=nanmax(FR_cond_sz8(:,id_cluster));
%             temp_max=max([FR_i_side1_lc1,FR_i_side2_lc1,FR_i_side1_lc2,FR_i_side2_lc2]);
            temp_max(temp_max==0)=1;
            BO_index_normbymax(id_cluster,id_ori,1)=(FR_i_side1_lc1-FR_i_side2_lc1)./temp_max;
            BO_index_normbymax(id_cluster,id_ori,2)=(FR_i_side1_lc2-FR_i_side2_lc2)./temp_max;
            BO_index_normbymax(id_cluster,id_ori,3)=(FR_i_side1_lc1+FR_i_side1_lc2-FR_i_side2_lc1-FR_i_side2_lc2)./(temp_max*2);
            LC_index_normbymax(id_cluster,id_ori,1)=(FR_i_side1_lc1-FR_i_side1_lc2)./temp_max;
            LC_index_normbymax(id_cluster,id_ori,2)=(FR_i_side2_lc1-FR_i_side2_lc2)./temp_max;
            LC_index_normbymax(id_cluster,id_ori,3)=(FR_i_side1_lc1+FR_i_side2_lc1-FR_i_side1_lc2-FR_i_side2_lc2)./(temp_max*2);

            BO_index_FR(id_cluster,id_ori,:)=[FR_i_side1_lc1,FR_i_side2_lc1,FR_i_side1_lc2,FR_i_side2_lc2];
            BO_index(id_cluster,id_ori,1)=(FR_i_side1_lc1-FR_i_side2_lc1)./(FR_i_side1_lc1+FR_i_side2_lc1);
            BO_index(id_cluster,id_ori,2)=(FR_i_side1_lc2-FR_i_side2_lc2)./(FR_i_side1_lc2+FR_i_side2_lc2);
            LC_index(id_cluster,id_ori,1)=(FR_i_side1_lc1-FR_i_side1_lc2)./(FR_i_side1_lc1+FR_i_side1_lc2);
            LC_index(id_cluster,id_ori,2)=(FR_i_side2_lc1-FR_i_side2_lc2)./(FR_i_side2_lc1+FR_i_side2_lc2);

            BO_index(id_cluster,id_ori,3)=(FR_i_side1_lc1+FR_i_side1_lc2-FR_i_side2_lc1-FR_i_side2_lc2)./(FR_i_side1_lc1+FR_i_side1_lc2+FR_i_side2_lc1+FR_i_side2_lc2);
            LC_index(id_cluster,id_ori,3)=(FR_i_side1_lc1+FR_i_side2_lc1-FR_i_side1_lc2-FR_i_side2_lc2)./(FR_i_side1_lc1+FR_i_side1_lc2+FR_i_side2_lc1+FR_i_side2_lc2);
        end
    end
end

p_anova_includedall=p_anova(Includedidx_all,:,:);
p_anova_side_includedall=p_anova_side(Includedidx_all,:,:);

BO_index_FR_included=BO_index_FR(Includedidx_all,:,:);
BO_index_normbymax(isnan(BO_index_normbymax))=0;
LC_index_normbymax(isnan(LC_index_normbymax))=0;
BO_index_included_normbymax=BO_index_normbymax(Includedidx_all,:,:);
LC_index_included_normbymax=LC_index_normbymax(Includedidx_all,:,:);
BO_index(isnan(BO_index))=0;
LC_index(isnan(LC_index))=0;
BO_index_included=BO_index(Includedidx_all,:,:);
LC_index_included=LC_index(Includedidx_all,:,:);
%%
% id_cond_reorder=[1,5,9,13,4,8,12,16,2,6,10,14,3,7,11,15];
% FR_cond_sz8_norm=FR_cond_sz8./repmat(max(FR_cond_sz8),16,1);
% FR_cond_sz8_reorder=FR_cond_sz8_norm(id_cond_reorder,:);
% FR_cond_ori2=FR_cond_sz8_norm(5:8,:);
% aa=FR_cond_sz8_reorder(1:8,Includedidx_all)'-FR_cond_sz8_reorder(9:16,Includedidx_all)';
% mean(aa)
% figure;
% cmap=brewermap(128,'*RdBu');
% for idx=1:3
%     subplot(1,3,idx)
%     if idx==1
%         image(FR_cond_sz8_reorder(1:8,Includedidx_all)','CDataMapping','scaled');
%         colormap(parula)
%     elseif idx==2
%         image(FR_cond_sz8_reorder(9:16,Includedidx_all)','CDataMapping','scaled');
%         colormap(parula)
%     elseif idx==3
%         image(FR_cond_ori2(1,Includedidx_all)'+FR_cond_ori2(3,Includedidx_all)'-FR_cond_ori2(2,Includedidx_all)'-FR_cond_ori2(4,Includedidx_all)','CDataMapping','scaled');        
% %         image(FR_cond_sz8_reorder(1:8,Includedidx_all)'-FR_cond_sz8_reorder(9:16,Includedidx_all)','CDataMapping','scaled');
%         colormap(cmap)
%     end
%     ax=gca;
%     if depth_isdeep==0
%         ax.YDir='Normal';
%     end
%     if idx<=2
%         ax.CLim=[0,1];  
%     else
%         ax.CLim=[-2,2];
%     end
%     ax.TickDir = 'out';
%     ax.Box='off';
%     ax.FontWeight='Bold';
%     ax.YAxis.LineWidth = 1;
%     ax.XAxis.LineWidth = 1;
%     xlabel('Orientation')
%     ylabel('Neuron ID (sorted by depth)')
%     title('Side')
% end
%% plot correlation of orientation tuning for size 4 and size 8
FR_ori_sz8_temp=reshape(FR_cond(2:2:32,:),4,4,N_cluster);
FR_ori_sz8=squeeze(nanmean(FR_ori_sz8_temp,1));
FR_ori_sz4_temp=reshape(FR_cond(1:2:31,:),4,4,N_cluster);
FR_ori_sz4=squeeze(nanmean(FR_ori_sz4_temp,1));

corr_ori=nan(N_cluster,1);
for i=1:N_cluster
    R=corrcoef(FR_ori_sz4(:,i),FR_ori_sz8(:,i));
    corr_ori(i)=R(1,2);    
end
% figure('Color',[1 1 1]);
[~,cluster.depthsorted_preori]=max(FR_ori_sz8,[],1);
FR_ori_norm_sz4=FR_ori_sz4./max(FR_ori_sz4,[],1);
FR_ori_norm_sz8=FR_ori_sz8./max(FR_ori_sz8,[],1);

cmap=brewermap(128,'Reds');
figure('Color',[1 1 1]);
% subplot(1,3,1)
% image(FR_ori_norm_sz4(:,Includedidx_all)','CDataMapping','scaled');
% colormap(cmap)
% hold on
% xLimits = get(gca,'XLim');
% if flag_exis_layer==1
%     for temp_i=1:length(cell_layer_idx_border_included)
%         plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
%         hold on
%     end 
% end
% ax=gca;
% if depth_isdeep==0
%     ax.YDir='Normal';
% end
% ax.CLim=[0,1];  
% ax.TickDir = 'out';
% ax.Box='off';
% ax.FontWeight='Bold';
% % ax.FontSize=sz_font;    
% ax.XTick=[1 2 3 4];
% ax.XTickLabel={'0','45','90','135'};
% % ax.YAxis.Visible='off';
% ax.YAxis.LineWidth = 1;
% ax.XAxis.LineWidth = 1;
% xlabel('Orientation')
% ylabel('Neuron ID (sorted by depth)')
% title('Size 4')
% 
% subplot(1,3,2)    
image(FR_ori_norm_sz8(:,Includedidx_all)','CDataMapping','scaled');
colormap(cmap)
hold on
xLimits = get(gca,'XLim');
if flag_exis_layer==1
    for temp_i=1:length(cell_layer_idx_border_included)
        plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
        hold on
    end 
end
ax=gca;
if depth_isdeep==0
    ax.YDir='Normal';
end
ax.CLim=[0,1];  
ax.TickDir = 'out';
ax.Box='off';
ax.FontWeight='Bold';
% ax.FontSize=sz_font;    
ax.XTick=[1 2 3 4];
ax.XTickLabel={'0','45','90','135'};
% ax.YAxis.Visible='off';
ax.YAxis.LineWidth = 1;
ax.XAxis.LineWidth = 1;
title('BO\_Size8')

% subplot(1,3,3)
% scatter(corr_ori(Includedidx_all),1:Ncell_included_all,'filled')
% hold on
% xLimits = get(gca,'XLim');
% if flag_exis_layer==1
%     for temp_i=1:length(cell_layer_idx_border_included)
%         plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
%         hold on
%     end 
% end
% title('cc (sz4 vs 8)')
% ax.TickDir = 'out';
% ax.Box='off';
% ax.FontWeight='Bold';
% ax.YAxis.LineWidth = 1;
% ax.XAxis.LineWidth = 1;
% xlabel('corr\_coeff')
% 
% xLimits = get(gca,'XLim');
% hold on
% if flag_exis_layer==1
%     for temp_i=1:length(cell_layer_idx_border_included)
%         plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
%         hold on
%     end 
% end
exportgraphics(gcf,[mysaveplotPath,'/','Orientation_BO','.pdf'],'ContentType','vector');  
%%
depthsort_included=round(cluster.depthsort(Includedidx_all'));
Nrows=ceil(max(depthsort_included)./500)*500;
FR_ori_sz8_norm_cmap=nan(Nrows,4);
FR_ori_sz4_norm_cmap=nan(Nrows,4);
FR_ori_norm_sz8_in=FR_ori_norm_sz8(:,Includedidx_all);
FR_ori_norm_sz4_in=FR_ori_norm_sz4(:,Includedidx_all);
for jjj=1:Ncell_included_all
    FR_ori_sz8_norm_cmap(depthsort_included(jjj)-4:depthsort_included(jjj)+4,:)=repmat(FR_ori_norm_sz8_in(:,jjj)',9,1);
    FR_ori_sz4_norm_cmap(depthsort_included(jjj)-4:depthsort_included(jjj)+4,:)=repmat(FR_ori_norm_sz4_in(:,jjj)',9,1);
end

figure('Color',[1 1 1]);
for id_sz=1:2
    subplot(1,2,id_sz)
    if id_sz==1
        image(1:4,1:Nrows,FR_ori_sz8_norm_cmap,'CDataMapping','scaled');
        colormap(cmap)
        title('BO\_Size8')
    else
        image(1:4,1:Nrows,FR_ori_sz4_norm_cmap,'CDataMapping','scaled');
        colormap(cmap)
        title('BO\_Size4')
    end
    hold on
    xLimits = get(gca,'XLim');
    if flag_exis_layer==1
        for id_layer=1:5
            plot(xLimits,[depth_mat(id_layer,2),depth_mat(id_layer,2)],'g--','LineWidth',layerline_width)
            hold on
        end 
    end
    ax=gca;
    if depth_isdeep==0
        ax.YDir='Normal';
    end
    ax.CLim=[0,1];  
    ax.TickDir = 'out';
    ax.Box='off';
    ax.FontWeight='Bold';
    % ax.FontSize=sz_font;   
    ax.YTick=0:100:Nrows;
    ax.XTick=[1 2 3 4];
    ax.XTickLabel={'0','45','90','135'};
    % ax.YAxis.Visible='off';
    ax.YAxis.LineWidth = 1;
    ax.XAxis.LineWidth = 1;
end

%%
cmap=brewermap(128,'Reds');
figure;
psth_cond_temp=squeeze(nanmean(psth_cond_all,1));
psth_cond_temp=psth_cond_temp./repmat(max(psth_cond_temp),N_Bins,1);
image(bincenter_psth,1:Ncell_included_all,psth_cond_temp(:,Includedidx_all)','CDataMapping','scaled');
hold on
colormap(cmap)
hold on
    for temp_i=1:length(cell_layer_idx_border_included)
        plot([0 1],[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
        hold on
    end 
ax=gca;
ax.CLim=[0,1.3];
ax.XLim=[0,0.3];
ax.YDir='Normal';
ax.TickDir = 'out';
ax.Box='off';
ax.FontWeight='Bold';
ax.FontSize=sz_font;    

% figure;
% for id_layer=1:4
%     psth_i=movmean(mean(psth_mean_temp(:,cell_layer_idx==id_layer&Includedidx'),2),100);
%     plot(bincenter_psth,psth_i-mean(psth_i(95:105)),'Color',colorlabel_layer(id_layer,:),'LineWidth',2);
%     hold on
% end
%% determine response latency
% data_s_included=data(:,1:500,Includedidx_all);
% latency=nan(Ncell_included_all,1);
% for i=1:Ncell_included_all
%     data_s_included_c=data_s_included(:,:,i);
%     isi_all=nan(size(data_s_included_c,1),1);
%     sptime=sp.st(sp.clu==cluster.depthsorted_id_includedall(i));      
%     mean_isi=nanmean(diff(sptime));   
%     burst=[];
%     for j=1:size(data_s_included_c,1)
%         spikestrain_c=find(data_s_included_c(j,:));
%         spikestrain_c=spikestrain_c./1000;
%         isi_j=nanmean(diff(spikestrain_c));
%         burst_temp=burstdetectps(spikestrain_c,mean_isi);
%         burst=[burst;burst_temp];        
%     end
%     burst_sig=burst(burst(:,6)>=8 & burst(:,1)>=0.05 & burst(:,1)<=0.3,:);
%     if size(burst_sig,1)>=3 
%         latency(i)=prctile(burst_sig(:,1),50);
%     end        
% end
% %%% plot latency
% figure;
% scatter(latency,1:Ncell_included_all,100,'k','filled') 
% xlim([0 0.3])
% hold on
% xLimits = get(gca,'XLim');
% if flag_exis_layer==1
%     for temp_i=1:length(cell_layer_idx_border_included)
%         plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
%         hold on
%     end 
% end
% ylim([0,Ncell_included_all])
%%
cluster.anovaa=p_anova_includedall;
cluster.cell_layer_idx_border_included=cell_layer_idx_border_included;
save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_cluster.mat']),'cluster','-v7.3');   

%% plot a hot map showing the location of significant side modulation neurons
% p_single_cmap=nan(8,Ncell_included);
% p_single_cmap([1,3,5,7],:)=squeeze(p_single(:,1,:));
% p_single_cmap([2,4,6,8],:)=squeeze(p_single(:,2,:));
%% calculate modulation index, which is the response ratio between (side1 - side2)/(side1 + side2)
%%% update: calculate based on (side1-side2)/max(FR)
% FR=nanmean(data_included,2)./psth_binSize;
% FR_cond=nan(32,Ncell_included);
% r_nc_all=nan(Ncell_included,Ncell_included,32);
% for idx_cond=1:32
%     FR_cond(idx_cond,:)=squeeze(nanmean(FR(id_condition==idx_cond,:,:),1))';
%     for i=1:Ncell_included-1
%         for j=i+1:Ncell_included
%             rtemp=corrcoef(FR(id_condition==idx_cond,:,i),FR(id_condition==idx_cond,:,j));
%             r_nc_all(i,j,idx_cond)=rtemp(1,2);
%             r_nc_all(j,i,idx_cond)=rtemp(1,2);
%         end
%     end
% end
% r_nc=nanmean(r_nc_all,3);
% FR_cond_temp=reshape(FR_cond,8,4,Ncell_included);
% FR_sz8=FR_cond_temp([2,4,8,6],:,:); %% 4cond * 4 ori* Ncell_included
% FR_sz4=FR_cond_temp([1,3,7,5],:,:); %% 4cond * 4 ori* Ncell_included
% BO_MI_lc1=nan(4,Ncell_included);
% BO_MI_lc2=nan(4,Ncell_included);
% FR_sz8_max=max(FR_sz8,[],[1,2]);
% for id_ori=1:4
%     BO_MI_lc1(id_ori,:)=(FR_sz8(1,id_ori,:)-FR_sz8(3,id_ori,:))./FR_sz8_max;
%     BO_MI_lc2(id_ori,:)=(FR_sz8(2,id_ori,:)-FR_sz8(4,id_ori,:))./FR_sz8_max;
% end
% for id_ori=1:4    
%     BO_MI_lc1(id_ori,:)=(FR_sz8(1,id_ori,:)-FR_sz8(3,id_ori,:))./(FR_sz8(3,id_ori,:)+std(FR(id_condition==(id_ori-1)*8+2 | id_condition==(id_ori-1)*8+8,1,:),'omitnan'));
%     BO_MI_lc2(id_ori,:)=(FR_sz8(2,id_ori,:)-FR_sz8(4,id_ori,:))./(FR_sz8(4,id_ori,:)+std(FR(id_condition==(id_ori-1)*8+4 | id_condition==(id_ori-1)*8+6,1,:),'omitnan'));
% end
% BO_MI_cmap=nan(8,Ncell_included);
% BO_MI_cmap([1,3,5,7],:)=BO_MI_lc1;
% BO_MI_cmap([2,4,6,8],:)=BO_MI_lc2;
% BO_MI_cmap(p_single_cmap>=0.01)=0;

% ax.Box='off';
% ax.XTickLabel={'Ori0\newlineLC1','\newline2','Ori45\newline LC1','\newline2','Ori90\newline LC1','\newline2','Ori135\newline LC1','\newline2'};
%% modulation ratio compariason between two lc conditions
% figure('Color',[1 1 1]);
% for id_ori=2
%     subplot(1,4,id_ori)
% %     scatterhist(BO_MI_lc1(id_ori,:),BO_MI_lc2(id_ori,:),'Group',cell_layer_idx','Kernel','on','Location','NorthEast','Direction','out','Color',colorlabel_layer,'LineWidth',2);
%     scatter(BO_MI_lc1(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig),BO_MI_lc2(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig),50,'k','filled');
%     hold on; axis square;xlim([-1,1]);ylim([-1,1])    
% end


%% side
% corr_corrected=corr_jk-corr_jk_shift;
% size(corr_corrected)
% figure('Position',[200,200,1000,1000],'Units', 'points');
% ax_pos_hot=[0.1, 0.1, 0.8, 0.8];
% ax=axes('Position', ax_pos_hot,'Units', 'points');
% cond_id=[2,4;8,6];
% modulation_corr=nan(4,Ncell_included);
% modulation_corr_1=nan(4,Ncell_included);
% modulation_corr_2=nan(4,Ncell_included);
% 
% modulation_corr_layer=nan(4,Ncell_included);
% modulation_corr_1_layer=nan(4,Ncell_included);
% modulation_corr_2_layer=nan(4,Ncell_included);
% 
% for id_ori=1:4
%     for id_side=1:3
%         if id_side<3
%             id_cond=cond_id(id_side,:)+(id_ori-1)*8;       
%             corr_condi=nanmax(nanmean(corr_corrected(:,:,:,id_cond),4),[],3);
%         else
%             corr_condi_1=nanmax(nanmean(corr_corrected(:,:,:,cond_id(1,:)+(id_ori-1)*8),4),[],3);
%             corr_condi_2=nanmax(nanmean(corr_corrected(:,:,:,cond_id(2,:)+(id_ori-1)*8),4),[],3);
%             corr_condi=(corr_condi_1-corr_condi_2)./(corr_condi_1+corr_condi_2);
%             modulation_corr(id_ori,:)=nanmean(corr_condi,1);
%             modulation_corr_1(id_ori,:)=nanmean(corr_condi_1,1);
%             modulation_corr_2(id_ori,:)=nanmean(corr_condi_2,1);
%             for id_layer=1:5
%                 modulation_corr_layer(id_ori,cell_layer_idx'==id_layer)=nanmean(corr_condi(cell_layer_idx'==id_layer,cell_layer_idx'==id_layer),1);
%                 modulation_corr_1_layer(id_ori,cell_layer_idx'==id_layer)=nanmean(corr_condi_1(cell_layer_idx'==id_layer,cell_layer_idx'==id_layer),1);
%                 modulation_corr_2_layer(id_ori,cell_layer_idx'==id_layer)=nanmean(corr_condi_2(cell_layer_idx'==id_layer,cell_layer_idx'==id_layer),1);
%             end
%         end
%         subplot(3,4,(id_side-1)*4+id_ori)
%         im1=image(corr_condi,'CDataMapping','scaled');
%         colormap(cmap)
%         if id_ori==1
%             colorbar
%         end
%         hold on
%         xLimits = get(gca,'XLim');
%         for temp_i=1:3
%             plot(xLimits,[cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],'k--','LineWidth',layerline_width)
%             hold on
%             plot([cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],xLimits,'k--','LineWidth',layerline_width)
%             hold on
%         end
%         ax=gca;
%         if depth_isdeep==0
%             ax.YDir='Normal';
%         end
%         if id_side<3
%             ax.CLim=[0,0.35]*0.001;
%         else
%             ax.CLim=[0,1000]*0.001;
%         end        
% %         if id_side<3
% %             ax.CLim=[0,0.15]*0.001;
% %         else
% %             ax.CLim=[0,0.06]*0.001;
% %         end       
%     end
% end
% filename=[recordingDate,'_',recordingSession,'_corr_side',filenameext];
% print('-painters','-dpdf', [mysaveplotPath,'\',filename,'.pdf'], '-bestfit');

%%
% %% local contrast
% corr_corrected=corr_jk-corr_jk_shift;
% size(corr_corrected)
% figure('Position',[200,200,1000,1000],'Units', 'points');
% ax_pos_hot=[0.1, 0.1, 0.8, 0.8];
% ax=axes('Position', ax_pos_hot,'Units', 'points');
% cond_id=[1,2,7,8;3,4,5,6];
% for id_ori=1:4
%     for id_lc=1:3
%         if id_lc<3
%             id_cond=cond_id(id_lc,:)+(id_ori-1)*8;   
%             corr_condi=nanmax(nanmean(corr_corrected(:,:,:,id_cond),4),[],3);
%         else
%             corr_condi=nanmax(nanmean(corr_corrected(:,:,:,cond_id(1,:)),4),[],3)-nanmax(nanmean(corr_corrected(:,:,:,cond_id(2,:)),4),[],3);
%         end                
%         
%         subplot(3,4,(id_lc-1)*4+id_ori)
%         im1=image(corr_condi,'CDataMapping','scaled');
%         colormap(cmap)
%         colorbar
%         hold on
%         xLimits = get(gca,'XLim');
%         for temp_i=1:3
%             plot(xLimits,[cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],'k--','LineWidth',layerline_width)
%             hold on
%             plot([cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],xLimits,'k--','LineWidth',layerline_width)
%             hold on
%         end
%         ax=gca;
%         ax.YDir='Normal';
%         ax.CLim=[0,0.15]*0.001;
%     end
% end
% 
% %% local size
% corr_corrected=corr_jk-corr_jk_shift;
% size(corr_corrected)
% figure('Position',[200,200,1000,1000],'Units', 'points');
% ax_pos_hot=[0.1, 0.1, 0.8, 0.8];
% ax=axes('Position', ax_pos_hot,'Units', 'points');
% cond_id=[1,3,5,7;2,4,6,8];
% for id_ori=1:4
%     for id_sz=1:3
%         if id_sz<3
%             id_cond=cond_id(id_sz,:)+(id_ori-1)*8;   
%             corr_condi=nanmax(nanmean(corr_corrected(:,:,:,id_cond),4),[],3);
%         else
%             corr_condi=nanmax(nanmean(corr_corrected(:,:,:,cond_id(1,:)),4),[],3)-nanmax(nanmean(corr_corrected(:,:,:,cond_id(2,:)),4),[],3);
%         end          
%         
%         subplot(3,4,(id_sz-1)*4+id_ori)
%         im1=image(corr_condi,'CDataMapping','scaled');
%         colormap(cmap)
%         colorbar
%         hold on
%         xLimits = get(gca,'XLim');
%         for temp_i=1:3
%             plot(xLimits,[cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],'k--','LineWidth',layerline_width)
%             hold on
%             plot([cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],xLimits,'k--','LineWidth',layerline_width)
%             hold on
%         end
%         ax=gca;
%         ax.YDir='Normal';
%         ax.CLim=[0,0.15]*0.001;
%     end
% end
% 
% %%
% corr_corrected=corr_jk-corr_jk_shift;
% size(corr_corrected)
% figure('Position',[200,200,1000,1000],'Units', 'points');
% ax_pos_hot=[0.1, 0.1, 0.8, 0.8];
% ax=axes('Position', ax_pos_hot,'Units', 'points');
% id_cond_pos=[1,2,3,4,7,8,5,6];
% for id_cond=1:8
%         corr_condi=nanmax(corr_corrected(:,:,:,id_cond),[],3);
%         corr_condi=[corr_condi;nan(1,size(corr_condi,2))];
%         subplot(2,4,id_cond_pos(id_cond))
%         im1=image(corr_condi,'CDataMapping','scaled');
%         colormap(cmap)
%         colorbar
%         hold on
%         xLimits = get(gca,'XLim');
%         for temp_i=1:3
%             plot(xLimits,[cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],'k--','LineWidth',layerline_width)
%             hold on
%             plot([cell_layer_idx_border(temp_i+1)+0.5,cell_layer_idx_border(temp_i+1)+0.5],xLimits,'k--','LineWidth',layerline_width)
%             hold on
%         end
%         ax=gca;
%         ax.YDir='Normal';
%         ax.CLim=[0,0.4]*0.001;
% end
% 


% %% plot histgram distribtuion of modulation index
% figure;
% histogram(Modulation_idx_maxori,10,'FaceColor',[0 0 0])
% ax=gca;
% ax.XLabel.String='Response ratio(R\_nonpref/R\_pref)';
% ax.YLabel.String='Number of Cells ';
% ax.XLabel.FontWeight='Bold';
% ax.YLabel.FontWeight='Bold';
% ax.Box='off';
% ax.Title.String=[recordingDate,'\_',recordingSession,', ','N=',num2str(Ncell_included)];
% filename=[recordingDate,'_',recordingSession,'_resoponse ratio'];
% print('-dpng', [mysaveplotPath,'\',filename,'.png'], '-r0');
% print('-painters','-dpdf', [mysaveplotPath,'\',filename,'.pdf'], '-r0'); 
% %% plot distribution of modulation index across depth

% figure('Color',[1 1 1],'WindowState', 'maximized');
% %%
% ax=gca;
% %% plot psth bos for each layer
% psth_bos=psth_prefside_maxori-psth_nprefside_maxori;
% figure('Color',[1 1 1],'WindowState', 'maximized');
% for id_layer=1:4
%     plot(bins,nanmean(psth_bos(cluster.depthsort>=depth_mat(id_layer,1)&cluster.depthsort<depth_mat(id_layer,2),:)),'Color',colorlabel2(id_layer,:),'Linewidth',2)
%     hold on
% end
% ax=gca;
% ax.XLabel.String='Time(ms)';
% ax.YLabel.String='Average responses';
% ax.XLabel.FontWeight='Bold';
% ax.YLabel.FontWeight='Bold';
% ax.Box='off';
% ax.TickDir = 'out';
% ax.PlotBoxAspectRatio = [1 1.5 1];
% ax.Title.String=[recordingDate,'\_',recordingSession,', ','\color[rgb]{0 0 0}N=',num2str(Ncell_included)];
% legend('5/6','4C','4B','2/3')
% legend('boxoff')
% filename=[recordingDate,'_',recordingSession,'_psth_layer'];
% print('-dpng', [mysaveplotPath,'\',filename,'.png'], '-r0');
% print('-painters','-dpdf', [mysaveplotPath,'\',filename,'.pdf'], '-r0');

% for i= 1:Ncell_included
%     i
%     idx=0;
%     sptime=sp.st(sp.clu==cluster.depthsorted_id(i));        
%     for ori_id=1:N_ori
%         for side_id=1:N_side
%             idx=idx+1;
%             event_on_time_currentcond=event_on_time(trialmat(:,c_ori)==ori_id&trialmat(:,c_sz)==1&trialmat(:,c_side)==side_id);
%             [psth, bins, rasterX, rasterY, spikeCounts, ba] = psthAndBA(sptime, event_on_time_currentcond, psth_window, psth_binSize);
%             % PSTH smoothing filter
%             gw = gausswin(round(smoothSize*6),3);
%             smWin = gw./sum(gw);        
%             % smooth ba
%             baSm = conv2(smWin,1,ba', 'same')'./psth_binSize; 
%             psthSm(idx,:,i) = mean(baSm);
%         end
%     end
% end
% psthSm_norm=psthSm./max(psthSm,[],[1 2]);
% psth_prefside_maxori=nan(Ncell_included,N_Bins);
% psth_nprefside_maxori=psth_prefside_maxori;
% for i= 1:Ncell_included
%     i
%     psth_prefside_maxori(i,:)=psthSm_norm((idx_maxori(i)-1)*2+idx_prefside(idx_maxori(i),i),:,i);
%     psth_nprefside_maxori(i,:)=psthSm_norm((idx_maxori(i)-1)*2+3-idx_prefside(idx_maxori(i),i),:,i);    
% end
% figure;
% plot(bins,nanmean(psth_prefside_maxori(p_maxori(:,1)<0.05,:)),'k','Linewidth',2)
% hold on
% plot(bins,nanmean(psth_nprefside_maxori(p_maxori(:,1)<0.05,:)),'k--','Linewidth',2)
% hold on
% ax=gca;
% ax.XLabel.String='Time(ms)';
% ax.YLabel.String='Average responses';
% ax.XLabel.FontWeight='Bold';
% ax.YLabel.FontWeight='Bold';
% ax.Box='off';
% legend('Pref Side','Non-pref Side')
% legend('boxoff')
% ax.Title.String=[recordingDate,'\_',recordingSession,', ','\color[rgb]{0 0 0}N=',num2str(Ncell_included),', sz4'];
% filename=[recordingDate,'_',recordingSession,'_psth_sz1'];
% print('-dpng', [myresultPath_psth,filename,'.png'], '-r300');





%% plot psth averaged across neurons, for side 1 and side 2
% figure;
% colorlabel2={'r','g','b','k'};
% for ori_id=1:N_ori
%     psth_bos=psthSm((ori_id-1)*2+1,:,:)-psthSm((ori_id-1)*2+2,:,:);
%     plot(bins,nanmean(psth_bos(1,:,:),3),colorlabel2{ori_id})
%     hold on
% end
% 
% for ori_id=1:N_ori
%     psth_bos=psthSm((ori_id-1)*2+1,:,:)-psthSm((ori_id-1)*2+2,:,:);
%     figure;
%     plot(bins,nanmean(psth_bos(1,:,depth_layer56_idx),3),'r')
%     hold on
%     plot(bins,nanmean(psth_bos(1,:,depth_layer4c_idx),3),'g')
%     hold on
%     plot(bins,nanmean(psth_bos(1,:,depth_layer4b_idx),3),'b')
%     hold on
%     plot(bins,nanmean(psth_bos(1,:,depth_layer23_idx),3),'m')
%     hold on
% end

%%


%         temp3=squeeze(psth_cond_included_sz8_norm(id_lc,id_ori,:,:)-psth_cond_included_sz8_norm(id_lc+2,id_ori,:,:));
%         im2=subplot(8,4,id_ori+8+(id_lc-1)*16);  
%         image(bincenter_psth(1:500),1:Ncell_included,temp3','CDataMapping','scaled');
%         colormap(im2,cmap1)    
%         ax=gca;
%         ax.YDir='Normal';
%         ax.CLim=[-0.8,0.8]; 
%         if id_ori>1
%             ax.YAxis.Visible='off';
%         end             
%         subplot(8,4,id_ori+12+(id_lc-1)*16);  
%         plot(bincenter_psth(1:500),psth_sz8(1,:),'b-','Linewidth',2)
%         hold on
%         plot(bincenter_psth(1:500),psth_sz8(2,:),'b:','Linewidth',2)
%         hold on        
%         plot(bincenter_psth(1:500),psth_sz8(1,:)-psth_sz8(2,:),'k-','Linewidth',3)
%         hold on
%         plot([bincenter_psth(1),bincenter_psth(500)],[0 0],'k:','Linewidth',1)
%         hold on        
%         ax=gca;
%         ax.XLim=[0,0.5];
%         ax.YLim=[-0.5,0.5];
%         if id_ori>1
%             ax.YAxis.Visible='off';
%         end      
