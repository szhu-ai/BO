maxlag=100;
T_tau=-maxlag:maxlag;
filenameext=['_',num2str(window_spikecount_cc(1)),'_',num2str(window_spikecount_cc(2))];    
load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_corr_jk',filenameext,'.mat']),'corr_jk');           
load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_corr_jk_jitter',filenameext,'.mat']),'corr_jk_jitter');

%% correlaton matrix is cell_j *cell_k* Ntau * Ncondition
% corr_corrected=corr_jk_jitter;
corr_corrected=corr_jk-corr_jk_jitter;
size(corr_corrected)
corr_all_mean=squeeze(nanmean(corr_corrected,4));
for i=1:size(corr_all_mean,1)
      corr_all_mean(i,i,:)=nan;
end
corr_all_flank=std(corr_all_mean(:,:,[1:50,end-49:end]),0,3);
corr_all_selected=corr_all_mean(:,:,find(T_tau==-10):find(T_tau==10));
[corr_all_peak,ccg_tau]=max(corr_all_selected,[],3);
ccg_sig_idx=corr_all_peak>=7*corr_all_flank;
ccg_sig_peak=corr_all_peak.*ccg_sig_idx;
ccg_sig_tau=ccg_tau.*ccg_sig_idx;
ccg_sig_tau(ccg_sig_tau==0)=nan;
ccg_sig_tau=ccg_sig_tau-11;
%%
cmap=brewermap(128,'*RdYlBu');
figure;
for id_ori=2
%     subplot(1,4,id_ori)
    im1=image(r_nc,'CDataMapping','scaled');
    colormap(cmap)
    colorbar
    ax=gca;
    if depth_isdeep==0
        ax.YDir='Normal';
    end
%     ax.CLim=[0,0.06]*0.001;
    hold on
    id_side_sig=find(p_all(:,1+(id_ori-1)*3)<0.01);
    if ~isempty(id_side_sig)
        plot([1,Ncell_included],[id_side_sig';id_side_sig'],'Color',[1 1 1],'LineWidth',1);
        hold on
    end
    if flag_exis_layer==1
        xLimits = get(gca,'XLim');
        for temp_i=1:length(cell_layer_idx_border)
            plot(xLimits,[cell_layer_idx_border(temp_i)+0.5,cell_layer_idx_border(temp_i)+0.5],'k--','LineWidth',layerline_width)
            hold on
            plot([cell_layer_idx_border(temp_i)+0.5,cell_layer_idx_border(temp_i)+0.5],xLimits,'k--','LineWidth',layerline_width)
            hold on
        end
    end
end

%%
cmap=brewermap(128,'*RdYlBu');
figure;
for id_ori=2
%     subplot(1,4,id_ori)
    im1=image(ccg_sig_peak,'CDataMapping','scaled');
    colormap(cmap)
    colorbar
    ax=gca;
    if depth_isdeep==0
        ax.YDir='Normal';
    end
%     ax.CLim=[0,0.06]*0.001;
    hold on
    id_side_sig=find(p_all(:,1+(id_ori-1)*3)<0.01);
    if ~isempty(id_side_sig)
        plot([1,Ncell_included],[id_side_sig';id_side_sig'],'Color',[1 1 1],'LineWidth',1);
        hold on
    end
    if flag_exis_layer==1
        xLimits = get(gca,'XLim');
        for temp_i=1:length(cell_layer_idx_border)
            plot(xLimits,[cell_layer_idx_border(temp_i)+0.5,cell_layer_idx_border(temp_i)+0.5],'k--','LineWidth',layerline_width)
            hold on
            plot([cell_layer_idx_border(temp_i)+0.5,cell_layer_idx_border(temp_i)+0.5],xLimits,'k--','LineWidth',layerline_width)
            hold on
        end
    end
end

%%
heccg_sig_peak_mean=nanmean(ccg_sig_peak,1);
ccg_sig_proportion=sum(ccg_sig_idx,1)./Ncell_included;
r_nc_mean=nanmean(r_nc);
figure('Color',[1 1 1]);
scatter(ccg_sig_peak_mean,1:Ncell_included,50,'k','filled');
figure('Color',[1 1 1]);
scatter(ccg_sig_proportion,1:Ncell_included,50,'k','filled');
figure('Color',[1 1 1]);
subplot(1,2,1)
scatter(BO_MI,ccg_sig_peak_mean,50,'k','filled');
subplot(1,2,2)
scatter(BO_MI_abs,ccg_sig_peak_mean,50,'k','filled');

figure('Color',[1 1 1]);
subplot(1,2,1)
scatter(BO_MI,ccg_sig_proportion,50,'k','filled');
subplot(1,2,2)
scatter(BO_MI_abs,ccg_sig_proportion,50,'k','filled');

figure('Color',[1 1 1]);
subplot(1,2,1)
scatter(BO_MI,r_nc_mean,50,'k','filled');
subplot(1,2,2)
scatter(BO_MI_abs,r_nc_mean,50,'k','filled');

figure('Color',[1 1 1]);
for id_ori=1:4
    for id_lc=1:2
        subplot(2,4,id_ori+(id_lc-1)*4)
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig,id_lc),ccg_sig_peak_mean(p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig),50,'k','filled');
        hold on 
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)>=p_sig,id_lc),ccg_sig_peak_mean(p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)>=p_sig),50,'r','filled');
        hold on 
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)<p_sig,id_lc),ccg_sig_peak_mean(p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)<p_sig),50,'g','filled');
        hold on    
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)<p_sig,id_lc),ccg_sig_peak_mean(p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)<p_sig),50,'b','filled');
        hold on  
        axis square
        xlim([-1,1])
    end
end
figure('Color',[1 1 1]);
for id_ori=1:4
    for id_lc=1:2
        subplot(2,4,id_ori+(id_lc-1)*4)
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig,id_lc),ccg_sig_proportion(p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig),50,'k','filled');
        hold on 
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)>=p_sig,id_lc),ccg_sig_proportion(p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)>=p_sig),50,'r','filled');
        hold on 
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)<p_sig,id_lc),ccg_sig_proportion(p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)<p_sig),50,'g','filled');
        hold on    
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)<p_sig,id_lc),ccg_sig_proportion(p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)<p_sig),50,'b','filled');
        hold on  
        axis square
        xlim([-1,1])
    end
end

figure('Color',[1 1 1]);
for id_ori=1:4
    for id_lc=1:2
        subplot(2,4,id_ori+(id_lc-1)*4)
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig,id_lc),r_nc_mean(p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)>=p_sig),50,'k','filled');
        hold on 
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)>=p_sig,id_lc),r_nc_mean(p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)>=p_sig),50,'r','filled');
        hold on 
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)<p_sig,id_lc),r_nc_mean(p_single(id_ori,1,:)>=p_sig & p_single(id_ori,2,:)<p_sig),50,'g','filled');
        hold on    
        scatter(BO_MI_lc_all(id_ori,p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)<p_sig,id_lc),r_nc_mean(p_single(id_ori,1,:)<p_sig & p_single(id_ori,2,:)<p_sig),50,'b','filled');
        hold on  
        axis square
        xlim([-1,1])
    end
end

%%
ccg_all_dist=nan(Ncell_included,Ncell_included);
for j=1:Ncell_included
    for k=1:Ncell_included
        ccg_all_dist(j,k)=abs(cluster.depthsort(j)-cluster.depthsort(k));
    end
end
ccg_sig_dist=ccg_all_dist.*ccg_sig_idx;
ccg_sig_dist_vec=reshape(ccg_sig_dist,[],1);
ccg_sig_tau_vec=reshape(ccg_sig_tau,[],1);
ccg_sig_peak_vec=reshape(ccg_sig_peak,[],1);

binranges=1:50:1000;
[bincounts,ind]=histc(ccg_sig_dist_vec,binranges);
ccg_dist_median_tau=nan(length(binranges),1);
ccg_dist_25_tau=nan(length(binranges),1);
ccg_dist_75_tau=nan(length(binranges),1);
ccg_dist_median_peak=nan(length(binranges),1);
ccg_dist_25_peak=nan(length(binranges),1);
ccg_dist_75_peak=nan(length(binranges),1);

for idx=1:length(binranges)
    ccg_dist_median_tau(idx)=nanmedian(abs(ccg_sig_tau_vec(ind==idx)));
    ccg_dist_25_tau(idx)=prctile(abs(ccg_sig_tau_vec(ind==idx)),25);
    ccg_dist_75_tau(idx)=prctile(abs(ccg_sig_tau_vec(ind==idx)),75);    
    ccg_dist_median_peak(idx)=nanmedian(abs(ccg_sig_peak_vec(ind==idx)));
    ccg_dist_25_peak(idx)=prctile(abs(ccg_sig_peak_vec(ind==idx)),25);
    ccg_dist_75_peak(idx)=prctile(abs(ccg_sig_peak_vec(ind==idx)),75);     
end
figure('Color',[1 1 1]);
scatter(binranges+25,ccg_dist_median_tau,80,'k','filled')
hold on
errorbar(binranges+25,ccg_dist_median_tau,ccg_dist_25_tau,ccg_dist_75_tau,'k','LineStyle','none','LineWidth',1)
ax=gca;
% ax.XLim=[0,360];
ax.Box='off';
ax.FontWeight='Bold';
ax.FontSize=16;
ax.TitleFontSizeMultiplier = 1;
ax.XLabel.String='vert. pair distance';  
ax.YLabel.String='avg. time lag';

figure('Color',[1 1 1]);
scatter(binranges+25,ccg_dist_median_peak,80,'k','filled')
hold on
errorbar(binranges+25,ccg_dist_median_peak,ccg_dist_25_peak,ccg_dist_75_peak,'k','LineStyle','none','LineWidth',1)
ax=gca;
% ax.XLim=[0,360];
ax.Box='off';
ax.FontWeight='Bold';
ax.FontSize=16;
ax.TitleFontSizeMultiplier = 1;
ax.XLabel.String='vert. pair distance';  
ax.YLabel.String='median peak';


%%    
figure;
for id_ori=1:4
    MI_c=BO_MI(id_ori,:);
    cross_MI=nan(Ncell_included,Ncell_included);
    for j=1:Ncell_included
        for k=1:Ncell_included
            cross_MI(j,k)=MI_c(j).*MI_c(k);
        end
    end
    cross_MI_sig=cross_MI.*ccg_sig_idx;
    subplot(1,4,id_ori)
    scatter(cross_MI_sig(:),ccg_sig_peak(:))
    hold on
end

%%

