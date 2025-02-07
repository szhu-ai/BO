%% plot psth averaged across neurons, for side 1 and side 2
psth_cond_all_reshaped=reshape(psth_cond_all,8,4,N_Bins,N_cluster);
% temp_max=max(psth_cond_all,[],[1,2],'omitnan');
% temp_max_sz8=mean(psth_cond_all(2:2:32,window_spikecount_idx,:),[1,2],'omitnan');
temp_max_sz8=max(psth_cond_all(2:2:32,window_spikecount_idx,:),[],[1,2],'omitnan');

idx_pref_ori_all=[0,0,0,0,1,2,0,0,0,1,0,0,4];
idx_pref_ori=idx_pref_ori_all(sessionidx);
psth_pref_ori=squeeze(psth_cond_all_reshaped(2:2:8,idx_pref_ori,:,:));
psth_side-
%%
psth_sz8_side_pref_allori=nan(2,N_Bins,N_cluster,4);
psth_sz8_lc_pref_allori=nan(2,N_Bins,N_cluster,4);
figure('Color',[1 1 1],'Position',[100 100 1400 1000]);   
for id_ori=1:4
    psth_sz8_side_pref_b=nan(2,N_Bins,N_cluster);
    psth_sz8_lc_pref_b=nan(2,N_Bins,N_cluster);
    psth_cond_sz8=squeeze(psth_cond_all_reshaped([2,4,6,8],id_ori,:,:));
    psth_cond_sz8_norm=psth_cond_sz8./repmat(temp_max_sz8,[4,N_Bins,1]);
    FR_i=squeeze(nanmean(psth_cond_sz8_norm(:,window_spikecount_idx,:),2));
    id_bo_sig=p_anova(:,id_ori,1)<0.05;
    id_lc_sig=p_anova(:,id_ori,2)<0.05;
    for id_cluster=1:N_cluster
        if Includedidx_all(id_cluster) && id_bo_sig(id_cluster)
            if FR_i(1,id_cluster)+FR_i(3,id_cluster)>=FR_i(2,id_cluster)+FR_i(4,id_cluster)
                psth_sz8_side_pref_b(:,:,id_cluster)=cat(1,mean(psth_cond_sz8_norm([1,3],:,id_cluster),1),mean(psth_cond_sz8_norm([2,4],:,id_cluster),1));
            elseif FR_i(1,id_cluster)+FR_i(3,id_cluster)<FR_i(2,id_cluster)+FR_i(4,id_cluster)
                psth_sz8_side_pref_b(:,:,id_cluster)=cat(1,mean(psth_cond_sz8_norm([2,4],:,id_cluster),1),mean(psth_cond_sz8_norm([1,3],:,id_cluster),1));
            end
        end
        if Includedidx_all(id_cluster) && id_lc_sig(id_cluster)    
            if FR_i(1,id_cluster)+FR_i(2,id_cluster)>=FR_i(3,id_cluster)+FR_i(4,id_cluster)
                psth_sz8_lc_pref_b(:,:,id_cluster)=cat(1,mean(psth_cond_sz8_norm([1,2],:,id_cluster),1),mean(psth_cond_sz8_norm([3,4],:,id_cluster),1));
            elseif FR_i(1,id_cluster)+FR_i(2,id_cluster)<FR_i(3,id_cluster)+FR_i(4,id_cluster)
                psth_sz8_lc_pref_b(:,:,id_cluster)=cat(1,mean(psth_cond_sz8_norm([3,4],:,id_cluster),1),mean(psth_cond_sz8_norm([1,2],:,id_cluster),1));
            end     
        end
    end
    for id_layer=N_layer-1:-1:1
        Includedidx_i=cell_layer_idx'==id_layer;
        N_ii=length(find(psth_sz8_side_pref_b(1,1,Includedidx_i)>=0));
        psth_cond_avg_b=squeeze(nanmean(psth_sz8_side_pref_b(:,:,Includedidx_i),3));
        psth_cond_std_b=nanstd(psth_sz8_side_pref_b(:,:,Includedidx_i),[],3)./sqrt(N_ii);
        psth_bo_avg_b=psth_cond_avg_b(1,:)-psth_cond_avg_b(2,:);
        psth_bo_std_b=nanstd(psth_sz8_side_pref_b(1,:,Includedidx_i)-psth_sz8_side_pref_b(2,:,Includedidx_i),[],3)./sqrt(N_ii);

        N_ii_lc=length(find(psth_sz8_lc_pref_b(1,1,Includedidx_i)>=0));
        psth_cond_lc_avg_b=squeeze(nanmean(psth_sz8_lc_pref_b(:,:,Includedidx_i),3));
        psth_cond_lc_std_b=nanstd(psth_sz8_lc_pref_b(:,:,Includedidx_i),[],3)./sqrt(N_ii_lc);
        psth_lc_avg_b=squeeze(nanmean(psth_sz8_lc_pref_b(1,:,Includedidx_i)-psth_sz8_lc_pref_b(2,:,Includedidx_i),3));
        psth_lc_std_b=nanstd(psth_sz8_lc_pref_b(1,:,Includedidx_i)-psth_sz8_lc_pref_b(2,:,Includedidx_i),[],3)./sqrt(N_ii_lc);

        subplot(2,4,id_ori)
        shadedErrorBar(bincenter_psth,psth_bo_avg_b,psth_bo_std_b,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
        hold on
        xlim([0,0.5])
        ylim([-0.1,0.5])
        subplot(2,4,id_ori+4)        
        shadedErrorBar(bincenter_psth,psth_lc_avg_b,psth_lc_std_b,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2,'LineStyle','-'}) 
        hold on
        xlim([0,0.5])
        ylim([-0.1,0.6])
    end

    psth_sz8_side_pref_allori(:,:,:,id_ori)=psth_sz8_side_pref_b;
    psth_sz8_lc_pref_allori(:,:,:,id_ori)=psth_sz8_lc_pref_b;
end
%%
figure('Color',[1 1 1],'Position',[100 100 1400 1000]);   
for id_layer=N_layer-1:-1:1
    Includedidx_i=cell_layer_idx'==id_layer;
    N_ii=length(find(psth_sz8_side_pref_allori(1,1,Includedidx_i,:)>=0));

    psth_cond_avg_b=squeeze(nanmean(psth_sz8_side_pref_allori(:,:,Includedidx_i,:),[3,4]));
    psth_cond_std_b=nanstd(psth_sz8_side_pref_allori(:,:,Includedidx_i,:),[],[3,4])./sqrt(N_ii);
    psth_bo_avg_b=psth_cond_avg_b(1,:)-psth_cond_avg_b(2,:);
    psth_bo_std_b=nanstd(psth_sz8_side_pref_allori(1,:,Includedidx_i,:)-psth_sz8_side_pref_allori(2,:,Includedidx_i,:),[],[3,4])./sqrt(N_ii);

    N_ii_lc=length(find(psth_sz8_lc_pref_allori(1,1,Includedidx_i,:)>=0));
    psth_cond_lc_avg_b=squeeze(nanmean(psth_sz8_lc_pref_allori(:,:,Includedidx_i,:),[3,4]));
    psth_cond_lc_std_b=nanstd(psth_sz8_lc_pref_allori(:,:,Includedidx_i,:),[],[3,4])./sqrt(N_ii_lc);
    psth_lc_avg_b=squeeze(nanmean(psth_sz8_lc_pref_allori(1,:,Includedidx_i,:)-psth_sz8_lc_pref_allori(2,:,Includedidx_i,:),[3,4]));
    psth_lc_std_b=nanstd(psth_sz8_lc_pref_allori(1,:,Includedidx_i,:)-psth_sz8_lc_pref_allori(2,:,Includedidx_i,:),[],[3,4])./sqrt(N_ii_lc);

    aa=[psth_cond_avg_s(1:2,:);psth_bo_avg_s(1,:);psth_cond_avg_s(3:4,:);psth_bo_avg_s(2,:);psth_cond_avg_b;psth_bo_avg_b;psth_cond_lc_avg_b;psth_lc_avg_b];
    ab=[psth_cond_std_s(1:2,:);psth_bo_std_s(1,:);psth_cond_std_s(3:4,:);psth_bo_std_s(2,:);psth_cond_std_b;psth_bo_std_b;psth_cond_lc_std_b;psth_lc_std_b];
    for idx=1:13
        if idx<=12
            subplot(5,3,idx)
            shadedErrorBar(bincenter_psth,aa(idx,:),ab(idx,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
        else
            subplot(5,3,12+ceil(id_layer./2))
            shadedErrorBar(bincenter_psth,psth_bo_avg_b,psth_bo_std_b,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
            hold on
            shadedErrorBar(bincenter_psth,psth_lc_avg_b,psth_lc_std_b,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2,'LineStyle','--'}) 
            hold on
        end

        hold on
        if mod(idx,3)==0 ||idx>12
            xlim([0,0.2])
        else
            xlim([-0.1,0.5])
        end
        ylim([-0.1,0.5])
        if id_layer==1 && idx==1
            xlabel('Time after Stim onset (s)')
            ylabel('BO\_LC1')
            title('Pref')
            legend({'L2/3','L4A/B','L4C','L5/6'},'Box','off');
        elseif idx==2
            title('Non\_pref')
        elseif idx==3
            title('Pref-non\_pref')
        elseif idx==4
            ylabel('BO\_LC2')
        elseif idx==7
            ylabel('BO\_LCall')
        elseif idx==10
            ylabel('LC')
        end
    end
end
filename='PSTH_BO_pref-nonpref';
exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')

%%
figure('Color',[1 1 1],'Position',[100 100 1400 500]);   
for id_ori=1:4
    psth_cond_sz8=squeeze(psth_cond_all_reshaped([2,4,6,8],id_ori,:,:));
    psth_cond_sz8_norm=movmean(psth_cond_sz8./repmat(temp_max_sz8,[4,N_Bins,1]),50,2);  

    for id_layer=N_layer-1:-1:1
%         Includedidx_i=Includedidx(:,id_ori) & cell_layer_idx'==id_layer;
        Includedidx_i=Includedidx_all & cell_layer_idx'==id_layer;
        N_ii=length(find(Includedidx_i));
        psth_cond_side1= mean(psth_cond_sz8_norm([1,3],:,Includedidx_i),1);
        psth_cond_side2= mean(psth_cond_sz8_norm([2,4],:,Includedidx_i),1);
        psth_cond_lc1= mean(psth_cond_sz8_norm([1,2],:,Includedidx_i),1);
        psth_cond_lc2= mean(psth_cond_sz8_norm([3,4],:,Includedidx_i),1);
        
        psth_bo=psth_cond_side1-psth_cond_side2;
        psth_bo_std=std(psth_bo,[],3)./sqrt(N_ii);
        psth_lc=psth_cond_lc1-psth_cond_lc2;
        psth_lc_std=std(psth_lc,[],3)./sqrt(N_ii);
        subplot(2,4,id_ori)
        shadedErrorBar(bincenter_psth,nanmean(psth_bo,3),psth_bo_std,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
        xlim([-0.1,0.5])
        ylim([-0.3,0.3])
        hold on
        if id_layer==1
            xlabel('Time after Stim onset (s)')
            ylabel('Side Modulation')
            title(['Ori\_',num2str(id_ori)])
        end
        subplot(2,4,id_ori+4)
        shadedErrorBar(bincenter_psth,nanmean(psth_lc,3),psth_lc_std,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
        hold on        
        xlim([-0.1,0.5])
        ylim([-0.3,0.3])

        if id_layer==1
            xlabel('Time after Stim onset (s)')
            ylabel('LC Modulation')
%             title(['Ori\_',num2str(id_ori)])
        end

    end
end
filename='PSTH_BO';
% exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')

%%
% figure('Color',[1 1 1],'Position',[100 100 1400 1000]);   
% for id_layer=1:N_layer-1
%     Includedidx_i=cell_layer_idx'==id_layer;
%     psth_cond_avg=squeeze(nanmean(psth_sz8_side_pref_s(:,:,Includedidx_i),3));
%     psth_cond_std=nanstd(psth_sz8_side_pref_s(:,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
%     psth_bo1_avg=squeeze(nanmean(psth_sz8_side_pref_s(1,:,Includedidx_i)-psth_sz8_side_pref_s(2,:,Includedidx_i),3));
%     psth_bo1_std=nanstd(psth_sz8_side_pref_s(1,:,Includedidx_i)-psth_sz8_side_pref_s(2,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
%     psth_bo2_avg=squeeze(nanmean(psth_sz8_side_pref_s(3,:,Includedidx_i)-psth_sz8_side_pref_s(4,:,Includedidx_i),3));
%     psth_bo2_std=nanstd(psth_sz8_side_pref_s(3,:,Includedidx_i)-psth_sz8_side_pref_s(4,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
% 
%     psth_cond_avg_b=squeeze(nanmean(psth_sz8_side_pref_b(:,:,Includedidx_i),3));
%     psth_cond_std_b=nanstd(psth_sz8_side_pref_b(:,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
%     psth_bo_avg_b=squeeze(nanmean(psth_sz8_side_pref_b(1,:,Includedidx_i)-psth_sz8_side_pref_b(2,:,Includedidx_i),3));
%     psth_bo_std_b=nanstd(psth_sz8_side_pref_b(1,:,Includedidx_i)-psth_sz8_side_pref_b(2,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
% 
%     subplot(3,7,id_layer)
%     shadedErrorBar(bincenter_psth,psth_cond_avg(1,:),psth_cond_std(1,:),'lineProps',{'Color',[1,0,0],'LineWidth',2})
%     hold on
%     shadedErrorBar(bincenter_psth,psth_cond_avg(2,:),psth_cond_std(2,:),'lineProps',{'Color',[0,0,1],'LineWidth',2})
%     hold on
%     title(['N=',num2str(length(find(Includedidx_i)))])
% %         ylim([0,0.8])
%     xlim([0,0.5])
%     subplot(3,7,5)
%     shadedErrorBar(bincenter_psth,psth_cond_avg(1,:),psth_cond_std(1,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
%     subplot(3,7,6)
%     shadedErrorBar(bincenter_psth,psth_cond_avg(2,:),psth_cond_std(2,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
%     subplot(3,7,7)
%     shadedErrorBar(bincenter_psth,psth_bo1_avg,psth_bo1_std,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
% 
%     subplot(3,7,id_layer+7)      
%     shadedErrorBar(bincenter_psth,psth_cond_avg(3,:),psth_cond_std(3,:),'lineProps',{'Color',[1,0,0],'LineWidth',2})
%     hold on
%     shadedErrorBar(bincenter_psth,psth_cond_avg(4,:),psth_cond_std(4,:),'lineProps',{'Color',[0,0,1],'LineWidth',2})
%     hold on
%     title(['N=',num2str(length(find(Includedidx_i)))])
% %         ylim([0,0.8])
%     xlim([0,0.5])
%     subplot(3,7,12)
%     shadedErrorBar(bincenter_psth,psth_cond_avg(3,:),psth_cond_std(3,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
%     subplot(3,7,13)
%     shadedErrorBar(bincenter_psth,psth_cond_avg(4,:),psth_cond_std(4,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
%     subplot(3,7,14)      
%     shadedErrorBar(bincenter_psth,psth_bo2_avg,psth_bo2_std,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
% %         ylim([0,0.8])
%     xlim([0,0.2])
% 
%     subplot(3,7,id_layer+14)      
%     shadedErrorBar(bincenter_psth,psth_cond_avg_b(1,:),psth_cond_std_b(1,:),'lineProps',{'Color',[1,0,0],'LineWidth',2})
%     hold on
%     shadedErrorBar(bincenter_psth,psth_cond_avg_b(2,:),psth_cond_std_b(2,:),'lineProps',{'Color',[0,0,1],'LineWidth',2})
%     hold on
%     title(['N=',num2str(length(find(Includedidx_i)))])
% %         ylim([0,0.8])
%     xlim([0,0.5])
%     subplot(3,7,19)
%     shadedErrorBar(bincenter_psth,psth_cond_avg_b(1,:),psth_cond_std_b(1,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
%     subplot(3,7,20)
%     shadedErrorBar(bincenter_psth,psth_cond_avg_b(2,:),psth_cond_std_b(2,:),'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
%     xlim([0,0.2])
%     subplot(3,7,21)      
%     shadedErrorBar(bincenter_psth,psth_bo_avg_b,psth_bo_std_b,'lineProps',{'Color',colorlabel_layer(id_layer,:),'LineWidth',2})
%     hold on
% %         ylim([0,0.8])
%     xlim([0.05,0.25])
% end


% 
%     figure('Color',[1 1 1]);
%     for id_layer=1:N_layer-1
%         Includedidx_i=Includedidx(:,id_ori)&cell_layer_idx'==id_layer;
%         psth_cond_avg=squeeze(nanmean(psth_cond_sz8_side_pref_s(:,:,Includedidx_i),3));
%         psth_cond_std=nanstd(psth_cond_sz8_side_pref_s(:,:,Includedidx(:,id_ori)&cell_layer_idx'==id_layer),[],3)./sqrt(length(find(Includedidx_i)));
%         subplot(3,4,id_layer)
%         shadedErrorBar(bincenter_psth,psth_cond_avg(1,:),psth_cond_std(1,:),'lineProps',{'Color',[1,0,0]})
%         hold on
%         shadedErrorBar(bincenter_psth,psth_cond_avg(2,:),psth_cond_std(2,:),'lineProps',{'Color',[0,0,1]})
%         hold on
%         title(['N=',num2str(length(find(Includedidx_i)))])
%         ylim([0,0.8])
%         xlim([0,0.5])
% 
%         subplot(3,4,id_layer+4)      
%         shadedErrorBar(bincenter_psth,psth_cond_avg(3,:),psth_cond_std(3,:),'lineProps',{'Color',[1,0,0]})
%         hold on
%         shadedErrorBar(bincenter_psth,psth_cond_avg(4,:),psth_cond_std(4,:),'lineProps',{'Color',[0,0,1]})
%         hold on        
%         ylim([0,0.8])
%         xlim([0,0.5])
% 
%         psth_cond_bo_pref=squeeze(nanmean(psth_bo_mod(:,:,Includedidx_i),3));
%         psth_bo_pref_std=nanstd(psth_bo_mod(:,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
%         subplot(3,4,id_layer+8)
%         shadedErrorBar(bincenter_psth,psth_cond_bo_pref,psth_bo_pref_std,'lineProps',{'Color',[0,0,0]})
%         hold on
%         plot([0,0.5],[0,0],'k--')
%         hold on
%         ylim([-0.1,0.3])
%         xlim([0,0.5])
%     end
%     filename=['BO index_timecourse_ori_',num2str(id_ori)];
%     exportgraphics(gcf,[mysaveplotPath,'/',filename,'.pdf'],'ContentType','vector'); 
% 
%     figure('Color',[1 1 1]);
%     for id_layer=1:N_layer-1
%         Includedidx_i=Includedidx(:,id_ori)&cell_layer_idx'==id_layer;
%         psth_cond_avg=squeeze(nanmean(psth_cond_sz8_lc_pref(:,:,Includedidx_i),3));
%         psth_cond_std=nanstd(psth_cond_sz8_lc_pref(:,:,Includedidx(:,id_ori)&cell_layer_idx'==id_layer),[],3)./sqrt(length(find(Includedidx_i)));
%         subplot(3,4,id_layer)
%         shadedErrorBar(bincenter_psth,psth_cond_avg(1,:),psth_cond_std(1,:),'lineProps',{'Color',[1,0,0]})
%         hold on
%         shadedErrorBar(bincenter_psth,psth_cond_avg(3,:),psth_cond_std(3,:),'lineProps',{'Color',[0,0,1]})
%         hold on
%         title(['N=',num2str(length(find(Includedidx_i)))])
%         ylim([0,0.8])
%         xlim([0,0.5])
% 
%         subplot(3,4,id_layer+4)      
%         shadedErrorBar(bincenter_psth,psth_cond_avg(2,:),psth_cond_std(2,:),'lineProps',{'Color',[1,0,0]})
%         hold on
%         shadedErrorBar(bincenter_psth,psth_cond_avg(4,:),psth_cond_std(4,:),'lineProps',{'Color',[0,0,1]})
%         hold on        
%         ylim([0,0.8])
%         xlim([0,0.5])
% 
%         psth_cond_lc_pref=squeeze(nanmean(psth_lc_mod(:,:,Includedidx_i),3));
%         psth_lc_pref_std=nanstd(psth_lc_mod(:,:,Includedidx_i),[],3)./sqrt(length(find(Includedidx_i)));
%         subplot(3,4,id_layer+8)
%         shadedErrorBar(bincenter_psth,psth_cond_lc_pref,psth_lc_pref_std,'lineProps',{'Color',[0,0,0]})
%         hold on
%         plot([0,0.5],[0,0],'k--')
%         hold on
%         ylim([-0.1,0.6])
%         xlim([0,0.5])
%     end
%     filename=['LC index_timecourse_ori_',num2str(id_ori)];
%     exportgraphics(gcf,[mysaveplotPath,'/',filename,'.pdf'],'ContentType','vector');  
%     close all

%%
cmap=parula(128); % 
% cmap=brewermap(128,'Reds');
temp_max=max(psth_cond_all(2:2:32,:,:),[],[1,2],'omitnan');
    figure('Color',[1 1 1],'Position',[100 100 1000 1000]);   
for id_ori=1:N_ori
    psth_cond_sz8=squeeze(psth_cond_all_reshaped([2,4,6,8],id_ori,:,:));
    psth_cond_sz8_norm=psth_cond_sz8./repmat(temp_max,[4,N_Bins,1]);
    for id_lc=1:2
        for id_side=1:2
            im1=nexttile(id_ori+(id_side-1)*4+(id_lc-1)*8);
            temp1=squeeze(psth_cond_sz8_norm(id_side+(id_lc-1)*2,:,Includedidx_all)); 
            image(bincenter_psth,1:Ncell_included_all,temp1','CDataMapping','scaled');
            colormap(im1,cmap)    
            ax=gca;
            if depth_isdeep==0
                ax.YDir='Normal';
            end
            ax.CLim=[0,1];
            if id_ori>1
                ax.YAxis.Visible='off';
            end     
            if id_ori==4 && id_side==2
                colorbar
            end 
            if id_side==1
                ax.XAxis.Visible='off';
            end 
            hold on
            xlim([-0.2,0.6])
            xLimits = get(gca,'XLim');
            if flag_exis_layer==1
                for temp_i=1:length(cell_layer_idx_border_included)
                    plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
                    hold on
                end 
            end

        end        
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exportgraphics(gcf,[mysaveplotPath,'/','psth.pdf'],'ContentType','vector');  
