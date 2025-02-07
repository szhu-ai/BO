close all
% method_label='_new';
method_label='';

%% load orientation ccg
path_ori='/Users/shushu/Dropbox/npix/code_new/CCG_V1/result';
path_bo_ori='/Users/shushu/Dropbox/npix/code_new/CCG_V1/BO_vs_ori';
load(fullfile(path_ori,['data_ori_',recordingDate,recordingSession,'_rmdouble_ccg_data_full_250_1000',method_label,'.mat']),'ccg_data_t_ori');
ccg_data_t_ori.ccg_control=fliplr(ccg_data_t_ori.ccg_control);
ccg_data_t_ori.ccg_norm=fliplr(ccg_data_t_ori.ccg_norm);
ccg_data_t_ori.ccg_norm_jitter=fliplr(ccg_data_t_ori.ccg_norm_jitter);
ccg_data_t_ori.geomean_nspike=sqrt(ccg_data_t_ori.pre_nspike.*ccg_data_t_ori.post_nspike);
ccg_data_t_ori.geomean_fr=sqrt(ccg_data_t_ori.pre_fr.*ccg_data_t_ori.post_fr);
ccg_data_t_ori.dist=sqrt((ccg_data_t_ori.xpos(:,2)-ccg_data_t_ori.xpos(:,1)).^2+(ccg_data_t_ori.ypos(:,2)-ccg_data_t_ori.ypos(:,1)).^2);

%% load BO ccg
flag.stim_size='all';
flag.start_time = 250; % ms, select spike train until 1000 ms after stimulus onset
flag.end_time = 1000; % ms, select spike train until 1000 ms after stimulus onset
if flag.stim_size==4
    sz_label='_sz4';
elseif flag.stim_size==8
    sz_label='_sz8';
else
    sz_label='_sz_all';   
end
load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_ccg_data_',num2str(flag.start_time),'_',num2str(flag.end_time),sz_label,method_label,'.mat']),'ccgs');   
load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_ccg_flag_',num2str(flag.start_time),'_',num2str(flag.end_time),sz_label,method_label,'.mat']),'flag');   

data_included_fr = squeeze(1000*mean(data(:, flag.start_time_index:flag.end_time_index, flag.Includedidx),[1,2]));
data_included_nspike = squeeze(sum(data(:, flag.start_time_index:flag.end_time_index, flag.Includedidx),[1,2]));
%% 
find(flag.Includedidx~=Includedidx_all)
ccg_mat=zeros(size(ccgs,1),flag.max_lag-flag.min_lag+1,'single');
ccg_norm_mat=zeros(size(ccgs,1),flag.max_lag-flag.min_lag+1,'single');
ccg_norm_jitter_mat=zeros(size(ccgs,1),flag.max_lag-flag.min_lag+1,'single');
ccg_data_t = struct;

for i=1:size(ccgs,1)
    if depth_isdeep==0
        ccg_norm_mat(i,:)=ccgs{i,1}.ccg_norm;
        ccg_norm_jitter_mat(i,:)=ccgs{i,1}.ccg_norm_jitter;
        pair_ids = flag.neuron_pair_all;   
    else
        ccg_norm_mat(i,:)=fliplr(ccgs{i,1}.ccg_norm);
        ccg_norm_jitter_mat(i,:)=fliplr(ccgs{i,1}.ccg_norm_jitter);  
        pair_ids = fliplr(flag.neuron_pair_all);    
    end
    ccg_mat(i,:)=ccg_norm_mat(i,:)-ccg_norm_jitter_mat(i,:);
end    
ccg_data_t.ccg_control = fliplr(ccg_mat);
ccg_data_t.ccg_norm=fliplr(ccg_norm_mat);
ccg_data_t.ccg_norm_jitter=fliplr(ccg_norm_jitter_mat);   
ccg_data_t.ccg_control_sm=ccg_data_t.ccg_control;
for i=1:size(ccg_data_t.ccg_control,1)
    ccg_data_t.ccg_control_sm(i,:)=conv(ccg_data_t.ccg_control(i,:),[0.05,0.25,0.4,0.25,0.05],'same');
end
[peaks, lpeak_lag] = max(ccg_data_t.ccg_control_sm,[],2);
[troughs, ltrough_lag] = min(ccg_data_t.ccg_control_sm,[],2);
ccg_data_t.config.flag = flag;    
ccg_data_t.peaks = peaks;
ccg_data_t.peaks_unsmooth=max(ccg_data_t.ccg_control,[],2);
ccg_data_t.troughs = troughs;        
ccg_data_t.peak_lag = single(lpeak_lag-101);
ccg_data_t.trough_lag = single(ltrough_lag-101);
%%
cluster.depthsorted_id_ccg=cluster.depthsorted_id(:,flag.Includedidx);
cluster.depthsorted_label_ccg=cluster.depthsorted_label(:,flag.Includedidx);
cluster.depthsort_ccg=cluster.depthsort(1,flag.Includedidx);
cluster.depthsorted_celltype_ccg=cluster.depthsorted_celltype(1,flag.Includedidx);
cluster.depthsorted_xpos_ccg=cluster.depthsorted_xpos(1,flag.Includedidx);
cluster.depthsorted_peakWF_ccg=cluster.depthsorted_peakWF(flag.Includedidx,:);

cluster.depthsorted_celllayer_ccg_ori=cell_layer_idx_ori(flag.Includedidx);
cluster.depthsorted_celllayer_full_ccg_ori=cell_layer_full_idx_ori(flag.Includedidx);
cluster.depthsorted_celllayer_ccg=cell_layer_idx(flag.Includedidx);
cluster.depthsorted_celllayer_full_ccg=cell_layer_full_idx(flag.Includedidx);

ccg_data_t.pre_id = pair_ids(:,1);
ccg_data_t.post_id = pair_ids(:,2);
ccg_data_t.pre_cluster_id = cluster.depthsorted_id_ccg(ccg_data_t.pre_id)';
ccg_data_t.post_cluster_id = cluster.depthsorted_id_ccg(ccg_data_t.post_id)';
ccg_data_t.pre_celltype = cluster.depthsorted_celltype_ccg(ccg_data_t.pre_id)';
ccg_data_t.post_celltype = cluster.depthsorted_celltype_ccg(ccg_data_t.post_id)';
ccg_data_t.pre_layerID = cluster.depthsorted_celllayer_ccg(ccg_data_t.pre_id)';
ccg_data_t.post_layerID = cluster.depthsorted_celllayer_ccg(ccg_data_t.post_id)';        
ccg_data_t.pre_layerID_full = cluster.depthsorted_celllayer_full_ccg(ccg_data_t.pre_id)';
ccg_data_t.post_layerID_full = cluster.depthsorted_celllayer_full_ccg(ccg_data_t.post_id)';        
ccg_data_t.pre_layerID_ori = cluster.depthsorted_celllayer_ccg_ori(ccg_data_t.pre_id)';
ccg_data_t.post_layerID_ori = cluster.depthsorted_celllayer_ccg_ori(ccg_data_t.post_id)';        
ccg_data_t.pre_layerID_full_ori = cluster.depthsorted_celllayer_full_ccg_ori(ccg_data_t.pre_id)';
ccg_data_t.post_layerID_full_ori = cluster.depthsorted_celllayer_full_ccg_ori(ccg_data_t.post_id)';        

ccg_data_t.pre_fr = data_included_fr(ccg_data_t.pre_id);
ccg_data_t.post_fr = data_included_fr(ccg_data_t.post_id);
ccg_data_t.geomean_fr=sqrt(ccg_data_t.pre_fr.*ccg_data_t.post_fr);
ccg_data_t.pre_nspike = data_included_nspike(ccg_data_t.pre_id);
ccg_data_t.post_nspike = data_included_nspike(ccg_data_t.post_id);
ccg_data_t.geomean_nspike=sqrt(ccg_data_t.pre_nspike.*ccg_data_t.post_nspike);

ccg_data_t.pre_BO_MI_LCall_oriall=BO_index_included(ccg_data_t.pre_id,:,3);
ccg_data_t.post_BO_MI_LCall_oriall=BO_index_included(ccg_data_t.post_id,:,3);
ccg_data_t.pre_BO_MI_LCall_oriall_normbymax=BO_index_included_normbymax(ccg_data_t.pre_id,:,3);
ccg_data_t.post_BO_MI_LCall_oriall_normbymax=BO_index_included_normbymax(ccg_data_t.post_id,:,3);

ccg_data_t.pre_BO_MI_LCall=mean(abs(BO_index_included(ccg_data_t.pre_id,:,3)),2);
ccg_data_t.post_BO_MI_LCall=mean(abs(BO_index_included(ccg_data_t.post_id,:,3)),2);
ccg_data_t.pre_BO_MI_LCall_normbymax=mean(abs(BO_index_included_normbymax(ccg_data_t.pre_id,:,3)),2);
ccg_data_t.post_BO_MI_LCall_normbymax=mean(abs(BO_index_included_normbymax(ccg_data_t.post_id,:,3)),2);

ccg_data_t.diff_BO_MI_LCall=mean(abs(BO_index_included(ccg_data_t.pre_id,:,3)-BO_index_included(ccg_data_t.post_id,:,3)),2);
ccg_data_t.diff_BO_MI_LCall_normbymax=mean(abs(BO_index_included_normbymax(ccg_data_t.pre_id,:,3)-BO_index_included_normbymax(ccg_data_t.post_id,:,3)),2);
ccg_data_t.geomean_BO_MI_LCall=sqrt(abs(ccg_data_t.pre_BO_MI_LCall).*abs(ccg_data_t.post_BO_MI_LCall));
ccg_data_t.geomean_BO_MI_LCall_normbymax=sqrt(abs(ccg_data_t.pre_BO_MI_LCall_normbymax).*abs(ccg_data_t.post_BO_MI_LCall_normbymax));

ccg_data_t.pre_LC_MI_sideall_oriall=LC_index_included(ccg_data_t.pre_id,:,3);
ccg_data_t.post_LC_MI_sideall_oriall=LC_index_included(ccg_data_t.post_id,:,3);
ccg_data_t.pre_LC_MI_sideall_oriall_normbymax=LC_index_included_normbymax(ccg_data_t.pre_id,:,3);
ccg_data_t.post_LC_MI_sideall_oriall_normbymax=LC_index_included_normbymax(ccg_data_t.post_id,:,3);

ccg_data_t.geomean_LC_MI_sideall=sqrt(mean(abs(ccg_data_t.pre_LC_MI_sideall_oriall),2).*mean(abs(ccg_data_t.post_LC_MI_sideall_oriall),2));
ccg_data_t.geomean_LC_MI_sideall_normbymax=sqrt(mean(abs(ccg_data_t.pre_LC_MI_sideall_oriall_normbymax),2).*mean(abs(ccg_data_t.post_LC_MI_sideall_oriall_normbymax),2));

bin_ccglag=-100:1:100;
bin_ccg_lead_id=find(bin_ccglag>=1&bin_ccglag<=13);
bin_ccg_lag_id=find(bin_ccglag>=-13&bin_ccglag<=-1);
integra_ccg_lead=sum(ccg_data_t.ccg_control(:,bin_ccg_lead_id),2);
integra_ccg_lag=sum(ccg_data_t.ccg_control(:,bin_ccg_lag_id),2);
ccg_data_t.CA=(integra_ccg_lead-integra_ccg_lag);


% ccg_data_t.r_sig=nan(1,length(ccg_data_t.pre_id));
% for i=1:length(ccg_data_t.pre_id)
%     temp=corrcoef(BO_index_FR_included(ccg_data_t.pre_id(i),id_max_ori,:),BO_index_FR_included(ccg_data_t.post_id(i),id_max_ori,:));
%     ccg_data_t.r_sig(i) = temp(1,2);
% end

pair_x_positions = cluster.depthsorted_xpos_ccg(pair_ids);    
pair_y_positions = cluster.depthsort_ccg(pair_ids);
ccg_data_t.xpos=pair_x_positions;
ccg_data_t.ypos=pair_y_positions;    
ccg_data_t.pair_distance  = sqrt(abs(pair_y_positions(:,1)-pair_y_positions(:,2)).^2+abs(pair_x_positions(:,1)-pair_x_positions(:,2)).^2);       
ccg_data_t.dist=sqrt((ccg_data_t.xpos(:,2)-ccg_data_t.xpos(:,1)).^2+(ccg_data_t.ypos(:,2)-ccg_data_t.ypos(:,1)).^2);
%%%%%% iterate over output struct and make everything verticle
fields = fieldnames(ccg_data_t);
for i = 1:length(fields)
    if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config')
        if size(ccg_data_t.(fields{i}), 1) == 1
            ccg_data_t.(fields{i}) = ccg_data_t.(fields{i})';
        end
        %ccg_data_t.(fields{i}) = single(ccg_data_t.(fields{i}));
    end
end

noise_distribution2 = [ccg_data_t.ccg_control_sm(:, 1:51), ccg_data_t.ccg_control_sm(:, 151:201)];
ccg_data_t.peaks_2nd=max(noise_distribution2,[],2);
ccg_data_t.noise_std2 = std(noise_distribution2,[],2, 'omitnan');
ccg_data_t.noise_mean2 = nanmean(noise_distribution2,2);
ccg_data_t.noise_median2 = prctile(abs(ccg_data_t.ccg_control),50,2);
ccg_data_t.peaks_ratio=(ccg_data_t.peaks-ccg_data_t.noise_mean2)./ccg_data_t.noise_std2;
ccg_data_t.sig_idx = (ccg_data_t.noise_std2>flag.sig_min_std ) &ccg_data_t.noise_median2>0&...
        (ccg_data_t.peaks>(flag.sig_num_stds*ccg_data_t.noise_std2 + ccg_data_t.noise_mean2)) & ...
        (abs(ccg_data_t.peak_lag) <= flag.sig_max_lag)&...
        ccg_data_t.geomean_fr>0.5&...
        (ccg_data_t.peaks_2nd<(3*ccg_data_t.noise_std2 + ccg_data_t.noise_mean2));
ccg_data = ccg_data_t;
fields = fieldnames(ccg_data);
for i = 1:length(fields)
    if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config')
        ccg_data.(fields{i}) = ccg_data.(fields{i})(ccg_data_t.sig_idx,:);
    end
end

ccg_data_t_s=ccg_data_t;
fields = fieldnames(ccg_data_t);
for i = 1:length(fields)
    if ~strcmp(fields{i},'cluster') && ~strcmp(fields{i},'config')
        ccg_data_t_s.(fields{i}) = ccg_data_t.(fields{i})(ccg_data_t.geomean_fr>0.5,:);
    end
end
disp("total number of pairs: " + length(ccg_data_t_s.sig_idx));
disp("number of significant pairs: " + sum(ccg_data_t.sig_idx));
save(fullfile(path_bo_ori,['flip_ccg_bo_full_',recordingDate,recordingSession,'_',num2str(flag.start_time),'_',num2str(flag.end_time),method_label,'.mat']),'ccg_data_t','-v7.3');
save(fullfile(path_bo_ori,['flip_ccg_ori_full_',recordingDate,recordingSession,'_',num2str(flag.start_time),'_',num2str(flag.end_time),method_label,'.mat']),'ccg_data_t_ori','-v7.3');

color_celltype=[0,0,0; 0,1,0;1 1 0; 1,0,0];
cmap=brewermap(256,'*RdBu'); %% prefer
layer_label={'5/6','4C','4A/B','2/3'};
layer_label_full={'6','5','4Cb','4Ca','4AB','3','2'};


%%
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(6,9,"TileSpacing","tight")
% for j=1:4
%     for k=(j+1):4
%         for id_pre_celltype=1:3
%             for id_post_celltype=1:3
%                 nexttile
%                 id_selected=ccg_data.pre_layerID==j&ccg_data.post_layerID==k&ccg_data.pre_celltype==id_pre_celltype&ccg_data.post_celltype==id_post_celltype;      
%                 id_selected_ori=ccg_data_ori.pre_layerID==j&ccg_data_ori.post_layerID==k&ccg_data_ori.pre_celltype==id_pre_celltype&ccg_data_ori.post_celltype==id_post_celltype;
%                 
%                 xx=ccg_data.peak_lag(id_selected);
%                 xx_ori=ccg_data_ori.peak_lag(id_selected_ori);
%        
%                 histogram(xx_ori,-11:2:11,'FaceColor',[0 0 0],'Normalization','count');
%                 hold on
%                 histogram(xx,-11:2:11,'FaceColor',color_celltype(id_pre_celltype+1,:),'Normalization','count');
%                 hold on
%                 yy=get(gca,'YLim');
%                 try
%                     p=ranksum(xx_ori,xx,'tail','both');
%                 catch
%                     p=nan;
%                 end
%                 if p<0.05
%                     title(['L',layer_label{j},'-L',layer_label{k},'\_',num2str(id_pre_celltype),'\_',num2str(id_post_celltype)],'FontSize',12,'Color',[0,0,1])
%                 else
%                     title(['L',layer_label{j},'-L',layer_label{k},'\_',num2str(id_pre_celltype),'\_',num2str(id_post_celltype)],'FontSize',12,'Color',[0,0,0])
%                 end
%                 hold on    
%                 if idx==1
%                     ylabel('Proportion of CCG')
%                     xlabel('Peak lag (ms)')
%                 end
%                 ax=gca;
%                 ax.Box='off';
%                 ax.YLim=[0,max([yy(2),10])];
%             end
%         end
%     end
% end
%%
% id_all=[ccg_data.pre_id;ccg_data.post_id];
% id_unique=unique(id_all);
% N_uniqueid=length(id_unique);
% aa=nan(N_uniqueid,1);
% for i=1:N_uniqueid
%     aa(i)=length(find(id_all==id_unique(i)));
% end
% 
% id_all=[ccg_data_ori.pre_id;ccg_data_ori.post_id];
% id_unique=unique(id_all);
% N_uniqueid=length(id_unique);
% ab=nan(N_uniqueid,1);
% for i=1:N_uniqueid
%     ab(i)=length(find(id_all==id_unique(i)));
% end
% 
% figure;
% histogram(aa)
% hold on
% histogram(ab)
% % N_lead_lag=nan(N_uniqueid,2);
% % for i=1:length(ccg_data.pre_id)
% %     N_lead_lag(i,1)=length(find(ccg_data.pre_id==i&ccg_data.peak_lag>0))+length(find(ccg_data.post_id==i&ccg_data.peak_lag<0));
% %     N_lead_lag(i,2)=length(find(ccg_data.pre_id==i&ccg_data.peak_lag<0))+length(find(ccg_data.post_id==i&ccg_data.peak_lag>0));
% % end
% 
% 
% %%
% [ccg_unique_id,~]=unique([ccg_data.pre_id;ccg_data.post_id]);
% figure;
% tiledlayout(5,5,"TileSpacing","tight")
% for idx=26:50
%     nexttile
%     id_selected=find(ccg_data.pre_id==ccg_unique_id(idx)|ccg_data.post_id==ccg_unique_id(idx));
%     for i=1:length(id_selected)
%         if ccg_data.peak_lag(id_selected(i))<0
%             plot([ccg_data.xpos(id_selected(i),1),ccg_data.xpos(id_selected(i),2)],[ccg_data.ypos(id_selected(i),1),ccg_data.ypos(id_selected(i),2)],'color',[0 0 1 0.5],'LineWidth',2)
%         elseif ccg_data.peak_lag(id_selected(i))==0
%             plot([ccg_data.xpos(id_selected(i),1),ccg_data.xpos(id_selected(i),2)],[ccg_data.ypos(id_selected(i),1),ccg_data.ypos(id_selected(i),2)],'color',[0 0 0 0.5],'LineWidth',2)
%         elseif ccg_data.peak_lag(id_selected(i))>0
%             plot([ccg_data.xpos(id_selected(i),1),ccg_data.xpos(id_selected(i),2)],[ccg_data.ypos(id_selected(i),1),ccg_data.ypos(id_selected(i),2)],'color',[1 0 0 0.5],'LineWidth',2)
%         end
%         hold on
%         scatter(ccg_data.xpos(id_selected(i),1),ccg_data.ypos(id_selected(i),1),'MarkerFaceColor',color_celltype(ccg_data.pre_celltype(i)+1,:));
%         hold on
%         scatter(ccg_data.xpos(id_selected(i),2),ccg_data.ypos(id_selected(i),2),'MarkerFaceColor',color_celltype(ccg_data.post_celltype(i)+1,:));
%         hold on
%     end
%     for id_layer=1:5
%         plot([10,70],[depth_layer_upper(id_layer),depth_layer_upper(id_layer)],'Color',[0 0 0 0.5],'LineWidth',2,'LineStyle','--')
%         hold on
%     end
%     ax=gca;
%     ax.XAxis.Visible='off';
%     ax.YAxis.Visible='off';
% 
% end
% 
% 
% %% plot ccg network
% figure;
% layer1=1;
% layer2=4;
% subplot(1,2,1)
% id_selected=find((ccg_data.pre_layerID==layer1&ccg_data.post_layerID==layer2)|(ccg_data.pre_layerID==layer2&ccg_data.post_layerID==layer1));
% for i=1:size(id_selected,1)
%     if ccg_data.peak_lag(id_selected(i))<0
%         plot([ccg_data.xpos(id_selected(i),1),ccg_data.xpos(id_selected(i),2)],[ccg_data.ypos(id_selected(i),1),ccg_data.ypos(id_selected(i),2)],'color',[0 0 1 0.3],'LineWidth',2)
%     elseif ccg_data.peak_lag(id_selected(i))==0 
%         plot([ccg_data.xpos(id_selected(i),1),ccg_data.xpos(id_selected(i),2)],[ccg_data.ypos(id_selected(i),1),ccg_data.ypos(id_selected(i),2)],'color',[0 0 0 0.3],'LineWidth',2)
%     elseif ccg_data.peak_lag(id_selected(i))>0 
%         plot([ccg_data.xpos(id_selected(i),1),ccg_data.xpos(id_selected(i),2)],[ccg_data.ypos(id_selected(i),1),ccg_data.ypos(id_selected(i),2)],'color',[1 0 0 0.3],'LineWidth',2)
%     end
%     hold on
% end
% for id_layer=1:5
%     plot([10,70],[depth_layer_upper(id_layer),depth_layer_upper(id_layer)],'Color',[0 0 0 0.5],'LineWidth',2,'LineStyle','--')
%     hold on
% end
% 
% subplot(1,2,2)
% id_selected=find((ccg_data_ori.pre_layerID==layer1&ccg_data_ori.post_layerID==layer2)|(ccg_data_ori.pre_layerID==layer2&ccg_data_ori.post_layerID==layer1));
% for i=1:size(id_selected,1)
%     if ccg_data_ori.peak_lag(id_selected(i))<0
%         plot([ccg_data_ori.xpos(id_selected(i),1),ccg_data_ori.xpos(id_selected(i),2)],[ccg_data_ori.ypos(id_selected(i),1),ccg_data_ori.ypos(id_selected(i),2)],'color',[0 0 1 0.3],'LineWidth',2)
%     elseif ccg_data_ori.peak_lag(id_selected(i))==0 
%         plot([ccg_data_ori.xpos(id_selected(i),1),ccg_data_ori.xpos(id_selected(i),2)],[ccg_data_ori.ypos(id_selected(i),1),ccg_data_ori.ypos(id_selected(i),2)],'color',[0 0 0 0.3],'LineWidth',2)
%     elseif ccg_data_ori.peak_lag(id_selected(i))>0 
%         plot([ccg_data_ori.xpos(id_selected(i),1),ccg_data_ori.xpos(id_selected(i),2)],[ccg_data_ori.ypos(id_selected(i),1),ccg_data_ori.ypos(id_selected(i),2)],'color',[1 0 0 0.3],'LineWidth',2)
%     end
%     hold on
% end
% for id_layer=1:5
%     plot([10,70],[depth_layer_upper(id_layer),depth_layer_upper(id_layer)],'Color',[0 0 0 0.5],'LineWidth',2,'LineStyle','--')
%     hold on
% end
% % [ccg_unique_id,idx]=unique([ccg_data_ori.pre_id;ccg_data_ori.post_id]);
% % unique_x_positions = cluster.depthsorted_xpos_includedall(ccg_unique_id);    
% % unique_y_positions = cluster.depthsort_includedall(ccg_unique_id);
% % unique_celltype = cluster.depthsorted_celltype_includedall(ccg_unique_id);
% % color_celltype=[0,0,0; 0,1,0;1 1 0; 1,0,0];
% % for i=1:length(ccg_unique_id)
% %     scatter(unique_x_positions(i),unique_y_positions(i),'MarkerFaceColor',color_celltype(unique_celltype(i)+1,:));
% % end
% %% plot single CCG example
% bin_ccglag=-100:1:100;
% sig_ccg_flip=ccg_data.ccg_control;
% sig_ccg_flip_sm=sig_ccg_flip;
% for i=1:size(sig_ccg_flip,1)
%     sig_ccg_flip_sm(i,:)=conv(sig_ccg_flip(i,:),[0.05,0.25,0.4,0.25,0.05],'same');
% end
% sig_ccg_flip_norm=sig_ccg_flip./repmat(max(sig_ccg_flip,[],2),1,size(sig_ccg_flip,2));
% sig_ccg_ori_flip=ccg_data_ori.ccg_norm-ccg_data_ori.ccg_norm_jitter;
% sig_ccg_norm=movmean(ccg_data.ccg_norm,2,2,"omitnan");
% sig_ccg_jitter=movmean(ccg_data.ccg_norm_jitter,2,2,"omitnan");
% 
% %%
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% for idx=61:80
%     subplot(4,5,idx-60)
%     plot(bin_ccglag,sig_ccg_ori_flip(idx,:),'k-')
%     hold on 
%     plot([0,0],get(gca,'YLim'),'k-')
%     xlim([-50,50])
%     title([num2str(ccg_data.pre_id(idx)),', ',num2str(ccg_data.post_id(idx))])
% end
% 
% 
% %% plot example CCG 
% id_selected_ori=find(ccg_data_ori.pre_layerID==3&ccg_data_ori.post_layerID==4);   
% id_selected=find(ccg_data.pre_layerID==3&ccg_data.post_layerID==4);            
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% for idx=1:min(20,length(id_selected_ori))
%     subplot(4,10,idx)
%     plot(bin_ccglag,sig_ccg_ori_flip(id_selected_ori(idx),:),'k-')
%     hold on 
%     plot([0,0],get(gca,'YLim'),'k-')
%     xlim([-20,20])
% end
% for idx=1:min(20,length(id_selected))
%     subplot(4,10,idx+20)
%     plot(bin_ccglag,sig_ccg_flip(id_selected(idx),:),'r-')
%     hold on 
%     plot([0,0],get(gca,'YLim'),'k-')
%     xlim([-20,20])
% end
% 
% %% plot CCG peak ratio
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(4,4,'TileSpacing','tight');
% sgtitle('Peaks\_sig','FontSize',20)
% for j=1:4
%     for k=1:4
%         if k>=j
%             nexttile((4-k)*4+j)
%             id_selected_all_ori=ccg_data_ori.pre_layerID==j&ccg_data_ori.post_layerID==k; 
%             xx1=ccg_data_ori.peaks(id_selected_all_ori);
%             histogram(xx1,0.001:0.002:0.04,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
% 
%             id_selected_all=ccg_data.pre_layerID==j&ccg_data.post_layerID==k;  
%             xx2=ccg_data.peaks(id_selected_all);
%             histogram(xx2,0.001:0.002:0.04,'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
%             try
%                 p=ranksum(xx1,xx2);
%             catch
%                 p=nan;
%             end
% %             plot([7,7],[0,0.2],'b-','LineWidth',1)
%             title(['L',layer_label{j},'-L',layer_label{k},',p=',num2str(p,2)],'FontSize',14)
%             ax=gca;
%             ax.Box='off';
%             if j>1
% %                 ax.YAxis.Visible='off';
%             end
%             ax.XLabel.String='Peaks';
%             if j~=k
% %                 ax.XAxis.Visible='off';
%             end
%             ax.YLim=[0,0.5];
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(4,4,'TileSpacing','tight');
% sgtitle('Peak\_ratio\_sig','FontSize',20)
% for j=1:4
%     for k=1:4
%         if k>=j
%             nexttile((4-k)*4+j)
%             id_selected_all_ori=ccg_data_ori.pre_layerID==j&ccg_data_ori.post_layerID==k;            
%             histogram(ccg_data_ori.peaks_ratio(id_selected_all_ori),1:1:20,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
% 
%             id_selected_all=ccg_data.pre_layerID==j&ccg_data.post_layerID==k;            
%             histogram(ccg_data.peaks_ratio(id_selected_all),1:1:20,'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
%             try
%                 p=ranksum(ccg_data_ori.peaks_ratio(id_selected_all_ori),ccg_data.peaks_ratio(id_selected_all));
%             catch
%                 p=nan;
%             end
%             plot([7,7],[0,0.2],'b-','LineWidth',1)
%             title(['L',layer_label{j},'-L',layer_label{k},',p=',num2str(p,2)],'FontSize',14)
%             ax=gca;
%             ax.Box='off';
%             if j>1
% %                 ax.YAxis.Visible='off';
%             end
%             ax.XLabel.String='Std above noise';
%             if j~=k
% %                 ax.XAxis.Visible='off';
%             end
%             ax.YLim=[0,0.5];
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(4,4,'TileSpacing','tight');
% sgtitle('Peak\_ratio\_all','FontSize',20)
% for j=1:4
%     for k=1:4
%         if k>=j
%             nexttile((4-k)*4+j)
%             id_selected_all_ori=ccg_data_t_ori.pre_layerID==j&ccg_data_t_ori.post_layerID==k;            
%             histogram(ccg_data_t_ori.peaks_ratio(id_selected_all_ori),1:1:20,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
% 
%             id_selected_all=ccg_data_t.pre_layerID==j&ccg_data_t.post_layerID==k;            
%             histogram(ccg_data_t.peaks_ratio(id_selected_all),1:1:20,'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
%             try
%                 p=ranksum(ccg_data_t_ori.peaks_ratio(id_selected_all_ori),ccg_data_t.peaks_ratio(id_selected_all));
%             catch
%                 p=nan;
%             end
%             plot([7,7],[0,0.2],'b-','LineWidth',1)
%             title(['L',layer_label{j},'-L',layer_label{k},',p=',num2str(p,2)],'FontSize',14)
%             ax=gca;
%             ax.Box='off';
%             if j>1
% %                 ax.YAxis.Visible='off';
%             end
%             ax.XLabel.String='Std above noise';
%             if j~=k
% %                 ax.XAxis.Visible='off';
%             end
%             ax.YLim=[0,0.5];
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(7,7,'TileSpacing','tight');
% for j=1:7
%     for k=1:7
%         if k>=j
%             nexttile((7-k)*7+j)
%             id_selected_all_ori=ccg_data_t_ori.pre_layerID_full==j&ccg_data_t_ori.post_layerID_full==k;            
%             histogram(ccg_data_t_ori.peaks_ratio(id_selected_all_ori),1:1:20,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
%             id_selected_all=ccg_data_t.pre_layerID_full==j&ccg_data_t.post_layerID_full==k;            
%             histogram(ccg_data_t.peaks_ratio(id_selected_all),1:1:20,'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
%             plot([7,7],[0,0.2],'r-','LineWidth',1)
%             title(['L',layer_label_full{j},'-L',layer_label_full{k}],'FontSize',16)
%             ax=gca;
%             ax.Box='off';
%             if j>1
% %                 ax.YAxis.Visible='off';
%             end
%             if j~=k
% %                 ax.XAxis.Visible='off';
%             end
%             ax.YLim=[0,0.5];
%         end
%     end
% end
% %% image plot of BO CCG, by layer
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(12,4,'TileSpacing','tight');
% for j=1:4
%     for k=1:4
%         if k>=j
%             id_selected=ccg_data.pre_layerID==j&ccg_data.post_layerID==k;
%             nexttile((4-k)*12+j)
%             histogram(ccg_data.peak_lag(id_selected),[-10.5:1:10.5]);
%             hold on
%             plot([0,0],get(gca,'YLim'),'k-','LineWidth',2)
%             hold on
%             title(['L',layer_label{j},'-L',layer_label{k}],'FontSize',16)
%             ax=gca;
%             ax.Box='off';
% %             ax.XAxis.Visible='off';
%             ax.XLim=[-20,20];
% 
%             nexttile((4-k)*12+j+4,[2,1])
%             image(bin_ccglag,1:sum(id_selected),sig_ccg_flip(id_selected,:),"CDataMapping","scaled")
%             colormap(jet)
%         %     colorbar
%             hold on
%             ax=gca;
%             ax.CLim=[-0.015,0.015];
%             ax.YDir='normal';
%             hold on    
%             plot([0,0],[0,sum(id_selected)],'k--','LineWidth',2)
%             hold on
%             yylim=get(gca,'YLim');
%             plot(bin_ccglag,yylim(2)/0.01*(0.002+mean(sig_ccg_flip(id_selected,:))),'k-','LineWidth',2);
% %             plot(bin_ccglag,yylim(2)/0.7*(0.15+mean(sig_ccg_flip_norm(id_selected,:))),'k-','LineWidth',3);
%             ax.Box='off';
% %             if j~=k
% %                 ax.XAxis.Visible='off';
% %                 %             ax.XTickLabel={};
% %             end
%             ax.YLim=[0,yylim(2)];
%             ax.XLim=[-50,50];
%         end
%     end
% end
% 
% %%% image plot of CCG_ori, by layer
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% tiledlayout(12,4,'TileSpacing','tight');
% for j=1:4
%     for k=1:4
%         if k>=j
%             id_selected=ccg_data_ori.pre_layerID==j&ccg_data_ori.post_layerID==k;
%             nexttile((4-k)*12+j)
%             histogram(ccg_data_ori.peak_lag(id_selected),[-10.5:1:10.5]);
%             hold on
%             plot([0,0],get(gca,'YLim'),'k-','LineWidth',2)
%             hold on
%             title(['L',layer_label{j},'-L',layer_label{k}],'FontSize',16)
%             ax=gca;
%             ax.Box='off';
% %             ax.XAxis.Visible='off';
%             ax.XLim=[-20,20];
% 
%             nexttile((4-k)*12+j+4,[2,1])
%             image(bin_ccglag,1:sum(id_selected),sig_ccg_ori_flip(id_selected,:),"CDataMapping","scaled")
%             colormap(jet)
%         %     colorbar
%             hold on
%             ax=gca;
%             ax.CLim=[-0.015,0.015];
%             ax.YDir='normal';
%             hold on    
%             plot([0,0],[0,sum(id_selected)],'k--','LineWidth',2)
%             hold on
%             yylim=get(gca,'YLim');
%             plot(bin_ccglag,yylim(2)/0.02*(0.005+mean(sig_ccg_ori_flip(id_selected,:))),'k-','LineWidth',2);
% %             plot(bin_ccglag,yylim(2)/0.7*(0.15+mean(sig_ccg_flip_norm(id_selected,:))),'k-','LineWidth',3);
%             ax.Box='off';
% %             if j~=k
% %                 ax.XAxis.Visible='off';
% %                 %             ax.XTickLabel={};
% %             end
%             ax.YLim=[0,yylim(2)];
%             ax.XLim=[-50,50];
%         end
%     end
% end
% %% plot mean CCG and peak lag distribution 
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% idx=0;
% for j=1:4
%     for k=j:4
%         idx=idx+1;
%         id_selected=ccg_data.pre_layerID==j&ccg_data.post_layerID==k;
%         id_selected_ori=ccg_data_ori.pre_layerID==j&ccg_data_ori.post_layerID==k;
%         subplot(6,5,idx)
%         xx1=ccg_data_ori.peak_lag(id_selected_ori);
%         xx2=ccg_data.peak_lag(id_selected);
%         histogram(xx1,-10.5:1:10.5,'FaceColor',[0 0 0],'Normalization','probability');
%         hold on
%         histogram(xx2,-10.5:1:10.5,'FaceColor',[1 0 0],'Normalization','probability');
%         hold on
%         try
%             p=ranksum(xx1,xx2,'tail','both');
%         catch
%             p=nan;
%         end
%         if p<0.05
%             title(['L',layer_label{j},'-L',layer_label{k},', p=',num2str(p,2)],'FontSize',12,'Color',[0,0,1])
%         else
%             title(['L',layer_label{j},'-L',layer_label{k},', p=',num2str(p,2)],'FontSize',12,'Color',[0,0,0])
%         end
%         hold on    
%         if idx==1
%             ylabel('Proportion of CCG')
%             xlabel('Peak lag (ms)')
%         end
%         ax=gca;
%         ax.Box='off';
% %         ax.YLim=[0,0.5];
%     
%         subplot(6,5,idx+10)
%         plot(bin_ccglag,nanmean(sig_ccg_flip(id_selected,:)),'Color',[1 0 0]);
%         hold on
%         plot(bin_ccglag,nanmean(sig_ccg_ori_flip(id_selected_ori,:)),'Color',[0 0 0]);
%         hold on
%         title(['L',layer_label{j},'-L',layer_label{k}])
%         hold on    
%         plot([0,0],get(gca,'YLim'),'k--')
%         xlim([-50,50])
%         ax=gca;
%         ax.Box='off';
% %         ax.YLim=[-0.005,0.015];
%         if idx==1
%             ylabel('Mean CCG (significant)')
%             xlabel('Peak lag (ms)')
%         end
% 
%         subplot(6,5,idx+20)
%         id_selected=ccg_data_t.pre_layerID==j&ccg_data_t.post_layerID==k;
%         id_selected_ori=ccg_data_t_ori.pre_layerID==j&ccg_data_t_ori.post_layerID==k;
%         plot(bin_ccglag,nanmean(ccg_data_t.ccg_control(id_selected,:)),'Color',[1 0 0]);
%         hold on
%         plot(bin_ccglag,nanmean(ccg_data_t_ori.ccg_control(id_selected_ori,:)),'Color',[0 0 0]);
%         hold on
%         title(['L',layer_label{j},'-L',layer_label{k}])
%         hold on    
%         plot([0,0],get(gca,'YLim'),'k--')
%         xlim([-50,50])
%         hold on
%         if idx==1
%             ylabel('Mean CCG (all)')
%             xlabel('Peak lag (ms)')
%         end
%         ax=gca;
%         ax.Box='off';
% %         ax.YLim=[-0.002,0.005];
%     end
% end
% %% plot full layer dependence
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]);
% for j=1:7
%     for k=1:7
%         if k>=j
%             subplot(7,7,(7-k)*7+j)
%             id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%             image(bin_ccglag,1:sum(id_selected),sig_ccg_flip(id_selected,:),"CDataMapping","scaled")
%             colormap(jet)
%     %         colorbar
%             xlim([-20,20])
%             hold on
%             ax=gca;
%             ax.CLim=[-0.015,0.015];
%             ax.YDir='normal';
%             title(['L',layer_label_full{j},'-L',layer_label_full{k}])
%             hold on    
%             plot([0,0],[1,sum(id_selected)],'k-','LineWidth',1)
%             hold on
%             yylim=get(gca,'YLim');
%             plot(bin_ccglag,yylim(2)/0.01*(0.002+mean(sig_ccg_flip(id_selected,:))),'k-','LineWidth',1);
% 
%             ax.Box='off';
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]); %%% mean ccg
% for j=1:7
%     for k=1:7
%         if k>=j
%             subplot(7,7,(7-k)*7+j)
%             id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%             plot(bin_ccglag,mean(sig_ccg_flip(id_selected,:)),'k','LineWidth',2);
%             xlim([-20,20])
%             hold on
%             ax=gca;
%             title(['L',layer_label_full{j},'-L',layer_label_full{k}],'FontSize',12)
%             hold on    
%             plot([-20,20],[0 0],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',2)
%             hold on
%             plot([0,0],[-0.002,0.01],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
%             hold on
%             ax.Box='off';
%             ax.YLim=[-0.002,0.01];
%             if j~=k
%                 ax.XAxis.Visible='off';
%                 %             ax.XTickLabel={};
%             end
%             if j>1
%               ax.YAxis.Visible='off';
%             end
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]); %%% peak lag distribution
% for j=1:7
%     for k=1:7
%         if k>=j
%             subplot(7,7,(7-k)*7+j)
%             id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%             histogram(ccg_data.peak_lag(id_selected),[-11:2:11]);
%             hold on
%             ax=gca;
%             title(['L',layer_label_full{j},'-L',layer_label_full{k}],'FontSize',12)
%             hold on    
%             yylim=get(ax,'YLim');
%             plot([0,0],[0,max([yylim(2),20])],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
%             hold on
%             ax.Box='off';
% %             ax.YLim=[-0.002,0.01];
%             ax.XLim=[-11,11];
%             if j~=k
% %                 ax.XAxis.Visible='off';
%                 ax.XTickLabel={};
%             end
%             if j>1
% %               ax.YAxis.Visible='off';
%             end
%         end
%     end
% end
% %%
% bin_ccg_lead_id=find(bin_ccglag>=0&bin_ccglag<=15);
% bin_ccg_lag_id=find(bin_ccglag>=-15&bin_ccglag<=0);
% integra_ccg_lead=sum(ccg_data.ccg_control(:,bin_ccg_lead_id),2);
% integra_ccg_lag=sum(ccg_data.ccg_control(:,bin_ccg_lag_id),2);
% CA=(integra_ccg_lead-integra_ccg_lag)./(abs(integra_ccg_lead)+abs(integra_ccg_lag));
% 
% integra_ccg_lead=sum(ccg_data_t.ccg_control(:,bin_ccg_lead_id),2);
% integra_ccg_lag=sum(ccg_data_t.ccg_control(:,bin_ccg_lag_id),2);
% CA_t=(integra_ccg_lead-integra_ccg_lag)./(abs(integra_ccg_lead)+abs(integra_ccg_lag));
% 
% integra_ccg_lead=sum(ccg_data_ori.ccg_control(:,bin_ccg_lead_id),2);
% integra_ccg_lag=sum(ccg_data_ori.ccg_control(:,bin_ccg_lag_id),2);
% CA_ori=(integra_ccg_lead-integra_ccg_lag)./(abs(integra_ccg_lead)+abs(integra_ccg_lag));
% 
% integra_ccg_lead=sum(ccg_data_t_ori.ccg_control(:,bin_ccg_lead_id),2);
% integra_ccg_lag=sum(ccg_data_t_ori.ccg_control(:,bin_ccg_lag_id),2);
% CA_t_ori=(integra_ccg_lead-integra_ccg_lag)./(abs(integra_ccg_lead)+abs(integra_ccg_lag));
% 
% if exist_sz==1
%     integra_ccg_lead=sum(aa.ccg_data.ccg_control(:,bin_ccg_lead_id),2);
%     integra_ccg_lag=sum(aa.ccg_data.ccg_control(:,bin_ccg_lag_id),2);
%     CA_sz4=(integra_ccg_lead-integra_ccg_lag)./(abs(integra_ccg_lead)+abs(integra_ccg_lag));
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]); 
% sgtitle('Correlogram Asymmetry','FontSize',16)
% for j=1:4
%     for k=1:4
%         if k>=j
%             subplot(4,4,(4-k)*4+j)
% 
%             id_selected_ori=ccg_data_ori.pre_layerID==j&ccg_data_ori.post_layerID==k;
%             histogram(CA_ori(id_selected_ori),[-1:0.2:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
% 
%             id_selected=ccg_data.pre_layerID==j&ccg_data.post_layerID==k;
%             histogram(CA(id_selected),[-1:0.2:1],'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
% 
%             ax=gca;
%             title(['L',layer_label{j},'-L',layer_label{k}],'FontSize',12)
%             hold on    
%             yylim=get(ax,'YLim');
%             plot([0,0],[0,1],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
%             hold on
%             ax.Box='off';
% %             ax.YLim=[-0.002,0.01];
%             ax.XLim=[-1,1];
%             ax.YLim=[0,0.8];
%             xlabel('CA')
%             if j==1
%                 ylabel('Proportion of CCGs')
%             end
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]); 
% sgtitle('Correlogram Asymmetry','FontSize',16)
% for j=1:4
%     for k=1:4
%         if k>=j
%             subplot(4,4,(4-k)*4+j)
%             id_selected_ori=ccg_data_t_ori.pre_layerID==j&ccg_data_t_ori.post_layerID==k;
%             histogram(CA_t_ori(id_selected_ori),[-1:0.2:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
% 
%             id_selected=ccg_data_t.pre_layerID==j&ccg_data_t.post_layerID==k;
%             histogram(CA_t(id_selected),[-1:0.2:1],'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
% 
%             ax=gca;
%             title(['L',layer_label{j},'-L',layer_label{k}],'FontSize',12)
%             hold on    
%             yylim=get(ax,'YLim');
%             plot([0,0],[0,1],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
%             hold on
%             ax.Box='off';
% %             ax.YLim=[-0.002,0.01];
%             ax.XLim=[-1,1];
%             ax.YLim=[0,0.8];
%             xlabel('CA')
%             if j==1
%                 ylabel('Proportion of CCGs')
%             end
%         end
%     end
% end
% 
% figure('Color',[1 1 1],'Position',[100,100,1400,1000]); 
% for j=1:7
%     for k=1:7
%         if k>=j
%             subplot(7,7,(7-k)*7+j)
% 
%             id_selected_ori=ccg_data_ori.pre_layerID_full==j&ccg_data_ori.post_layerID_full==k;
%             histogram(CA_ori(id_selected_ori),[-1:0.2:1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on
% 
%             id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%             histogram(CA(id_selected),[-1:0.2:1],'FaceColor',[1 0 0],'Normalization','probability');
%             hold on
% 
%             ax=gca;
%             title(['L',layer_label_full{j},'-L',layer_label_full{k}],'FontSize',12)
%             hold on    
%             yylim=get(ax,'YLim');
%             plot([0,0],[0,max([yylim(2),1])],'Color',[0.5 0.5 0.5],'LineStyle','--','LineWidth',2)
%             hold on
%             ax.Box='off';
% %             ax.YLim=[-0.002,0.01];
%             ax.XLim=[-1,1];
%             if j~=k
% %                 ax.XAxis.Visible='off';
%                 ax.XTickLabel={};
%             end
%             if j>1
% %               ax.YAxis.Visible='off';
%             end
%         end
%     end
% end
% %%
% for id_pre_type=[1,3]
%     for id_post_type=[1,3]
%         figure('Position',[100,100,800,500]);
%         for j=1:7
%            for k=1:7
%                 subplot(7,7,(7-j)*7+k)
%                 id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k&ccg_data.pre_celltype==id_pre_type&ccg_data.post_celltype==id_post_type;
%                 histogram(CA(id_selected),[-1.1:0.2:1.1],'FaceColor',[0.5 0.5 0.5]);
%                 hold on
%                 plot([0,0],[0,5],'k--','LineWidth',2)
%                 hold on
%                 ax=gca;
%                 ax.Box='off';
%                 title(['L',layer_label{j},'-L',layer_label{k}])
%                 hold on    
%            end
%         end
%         sgtitle(['Pretype=',num2str(id_pre_type),', Posttype=',num2str(id_post_type)])
%     end
% end
% %%
% if exist_sz==1
%     figure;
%     subplot(2,1,1)
%     histogram(aa.ccg_data.peaks,0:0.002:0.06,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%     hold on
%     histogram(ccg_data.peaks,0:0.002:0.06,'FaceColor',[0.5 0 0],'Normalization','probability');
%     hold on
%     ranksum(ccg_data.peaks,aa.ccg_data.peaks)
%     subplot(2,1,2)
%     histogram(aa.ccg_data.peak_lag,-10:1:10,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%     hold on
%     histogram(ccg_data.peak_lag,-10:1:10,'FaceColor',[0.5 0 0],'Normalization','probability');
%     hold on
%     ranksum(ccg_data.peak_lag,aa.ccg_data.peak_lag)
% 
%     %%%%%%%% 
%     figure;
%     for j=1:7
%        for k=1:7
%             subplot(7,7,(7-j)*7+k)
%             id_selected=aa.ccg_data.pre_layerID_full==j&aa.ccg_data.post_layerID_full==k;
%             histogram(aa.ccg_data.peaks(id_selected),0:0.002:0.06,'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%             hold on  
%     
%             id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%             histogram(ccg_data.peaks(id_selected),0:0.002:0.06,'FaceColor',[0.5 0 0],'Normalization','probability');
%             hold on
%           
%             hold on
%             ax=gca;
%             ax.Box='off';
%             title(['layer',layer_label{j},'-layer',layer_label{k}])
%             hold on    
%        end
%     end
% 
%     %%%%%%%%%%%%%
%     % figure;
%     % for j=1:7
%     %    for k=1:7
%     %         subplot(7,7,(7-j)*7+k)
%     %         id_selected=aa.ccg_data.pre_layerID_full==j&aa.ccg_data.post_layerID_full==k;
%     %         histogram(CA_sz4(id_selected),[-1.1:0.2:1.1],'FaceColor',[0.5 0.5 0.5],'Normalization','probability');
%     %         hold on  
%     % 
%     %         id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%     %         histogram(CA(id_selected),[-1.1:0.2:1.1],'FaceColor',[0.5 0 0],'Normalization','probability');
%     %         hold on
%     %       
%     %         hold on
%     %         ax=gca;
%     %         ax.Box='off';
%     %         title(['layer',layer_label{j},'-layer',layer_label{k}])
%     %         hold on    
%     %    end
%     % end
%     figure;
%     for j=1:7
%        for k=1:7
%             subplot(7,7,(7-j)*7+k)
%             id_selected=aa.ccg_data.pre_layerID_full==j&aa.ccg_data.post_layerID_full==k;
%             histogram(CA_sz4(id_selected),[-1.1:0.2:1.1],'FaceColor',[0.5 0.5 0.5]);
%             hold on  
%     
%             id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%             histogram(CA(id_selected),[-1.1:0.2:1.1],'FaceColor',[0.5 0 0]);
%             hold on
%           
%             hold on
%             ax=gca;
%             ax.Box='off';
%             title(['layer',layer_label{j},'-layer',layer_label{k}])
%             hold on    
%        end
%     end
% end
% %%
% figure;
% histogram(ccg_data_t.r_sig,-1:0.1:1,'Normalization','probability','FaceColor',[0.5,0.5,0.5])
% hold on
% histogram(ccg_data.r_sig,-1:0.1:1,'Normalization','probability','FaceColor',[0.5,0,0])
% hold on
% 
% % % scatter(ccg_data.r_sig,ccg_data.peaks,10,'filled');
% 
% %%
% peaks_nonsig=ccg_data_t.peaks;
% peaks_nonsig(peaks_nonsig>=0.06)=0.06;
% peaks_sig=ccg_data.peaks;
% peaks_sig(peaks_sig>=0.06)=0.06;
% figure;
% scatter(ccg_data_t.r_sig,peaks_nonsig,10,'filled','MarkerFaceColor',[0.5 0.5 0.5]);
% hold on
% scatter(ccg_data.r_sig,peaks_sig,10,'filled','MarkerFaceColor',[1,0,0]);
% 
% %%
% figure;
% for j=1:7
%    for k=1:7
%         subplot(7,7,(7-j)*7+k)
% %         id_selected=ccg_data_t.pre_layerID_full==j&ccg_data_t.post_layerID_full==k;
% %         scatter(ccg_data_t.pre_ModIdx_2(id_selected),ccg_data_t.post_ModIdx_2(id_selected),10,'filled');
% %         hold on
% 
%         id_selected=ccg_data_t.pre_layerID_full==j&ccg_data_t.post_layerID_full==k;
%         scatter(ccg_data_t.r_sig(id_selected),peaks_nonsig(id_selected),10,'filled','MarkerFaceColor',[0.5 0.5 0.5]);
%         hold on
%         id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
%         scatter(ccg_data.r_sig(id_selected),peaks_sig(id_selected),30,'filled','MarkerFaceColor',[1 0 0]);
%         ax=gca;
%         xlim([-1,1]);
%         ylim([0,0.06]);
%         title(['layer',layer_label{j},'-layer',layer_label{k}])
%         hold on    
%    end
% end
% %%
% % figure;
% % for j=1:7
% %    for k=1:7
% %         subplot(7,7,(7-j)*7+k)
% %         id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
% %         scatter(abs(ccg_data.pre_ModIdx_2(id_selected)).*abs(ccg_data.post_ModIdx_2(id_selected)),ccg_data.peaks(id_selected),20,'filled','MarkerFaceColor',[0.5 0.5 0.5]);
% %         hold on
% %         scatter(abs(ccg_data.pre_ModIdx_1(id_selected)).*abs(ccg_data.post_ModIdx_1(id_selected)),ccg_data.peaks(id_selected),20,'filled','MarkerFaceColor',[1 0 0]);
% %         ax=gca;
% % %         xlim([-1,1]);
% % %         ylim([0,0.06]);
% %         title(['layer',layer_label{j},'-layer',layer_label{k}])
% %         hold on    
% %    end
% % end
% %%
% % figure;
% % for j=1:7
% %    for k=1:7
% %         subplot(7,7,(7-j)*7+k)
% %         id_selected=ccg_data_t.pre_layerID_full==j&ccg_data_t.post_layerID_full==k;
% % %         histogram(ccg_data_t.pre_ModIdx_3(id_selected),[-1.1:0.2:1.1],'FaceColor',[0,0,0],'Normalization','probability');
% % %         hold on
% %         histogram(ccg_data_t.post_ModIdx_1(id_selected),[-1.1:0.2:1.1],'FaceColor',[0.5,0.5,0.5],'Normalization','probability');
% %         hold on
% % 
% %         id_selected=ccg_data.pre_layerID_full==j&ccg_data.post_layerID_full==k;
% % %         histogram(ccg_data.pre_ModIdx_3(id_selected),[-1.1:0.2:1.1],'FaceColor',[1,0,0],'Normalization','probability');
% % %         hold on
% %         histogram(ccg_data.post_ModIdx_1(id_selected),[-1.1:0.2:1.1],'FaceColor',[0,0,1],'Normalization','probability');
% %         hold on
% %         ax=gca;
% %         title(['layer',layer_label{j},'-layer',layer_label{k}])
% %         hold on    
% %    end
% % end