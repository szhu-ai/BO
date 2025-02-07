%% plot psth averaged across neurons, for side 1 and side 2
close all
psth_cond_included=psth_cond_all(:,:,Includedidx_all);

Ncell_included_all=sum(Includedidx_all);
psth_cond_included_temp=reshape(psth_cond_included,8,4,1300,Ncell_included_all);
psth_cond_included_sz8=psth_cond_included_temp([2,4,6,8],:,:,:);
psth_cond_included_sz4=psth_cond_included_temp([1,3,5,7],:,:,:);
psth_cond_included_sz_all=cat(5,psth_cond_included_sz8,psth_cond_included_sz4);

FR_cond_reshape=reshape(FR_cond,8,4,N_cluster);
FR_sz8=FR_cond_reshape([2,4,6,8],:,Includedidx_all);
FR_sz4=FR_cond_reshape([1,3,5,7],:,Includedidx_all);
FR_sz_all=cat(4,FR_sz8,FR_sz4);
FR_mat_sz8=FR_mat_reshape_sz8(:,:,:,Includedidx_all);
FR_mat_sz4=FR_mat_reshape_sz4(:,:,:,Includedidx_all);
FR_mat_sz_all=cat(5,FR_mat_sz8,FR_mat_sz4);
delay_sys=0.03;
colorscheme1=brewermap(8,'Dark2');
%% plot psth for each single neuron
square_ori=[0,45,90,135];
colorlabel_side=[colorscheme1(4,:); colorscheme1(5,:)];
unitlabel={'MU','SU'};
% colorlabel_side=[[0 0 0]; [0 0 0]];
% delay_sys=0;
for id_sz=1
    for id_ori=1
        figure('Color',[1 1 1],'Position',[100 100 600 600],'Visible','on');     
        for id_lc=1:2
            %% raster plot
            for i= 1:Ncell_included_all %, example[98,96,86,60]
                sptime=sp.st(sp.clu==cluster.depthsorted_id_includedall(i)); 
                event_on_time_currentcond_side1=event_on_time(trialmat(:,c_ori)==id_ori&trialmat(:,c_side)==1&trialmat(:,c_lc)==id_lc&trialmat(:,c_sz)==3-id_sz); 
                [~, ~, rasterx_side1, rasterY_side1, ~, ~] = psthAndBA(sptime, event_on_time_currentcond_side1, [0,1], psth_binSize);
                event_on_time_currentcond_side2=event_on_time(trialmat(:,c_ori)==id_ori&trialmat(:,c_side)==2&trialmat(:,c_lc)==id_lc&trialmat(:,c_sz)==3-id_sz); 
                [~, ~, rasterx_side2, rasterY_side2, ~, ~] = psthAndBA(sptime, event_on_time_currentcond_side2, [0,1], psth_binSize);
                ax=subplot(1,4,(id_lc-1)*2+1);                 
                plot(rasterx_side1-delay_sys,rasterY_side1+N_repetition*(i-1),'Color',colorlabel_side(1,:),'LineWidth',1);
                hold on
%                 ax.YLim=[800,2400];
                ax.XLim=[-0.2,0.5];
                ax.TickDir='out';
                ax.Box='off';
                ax.TickLength=[0.02,0.02];
                ax.FontSize=6;
%                 ax.YTickLabel={};
                ax=subplot(1,4,(id_lc-1)*2+2);                 
                plot(rasterx_side2-delay_sys,rasterY_side2+N_repetition*(i-1),'Color',colorlabel_side(2,:),'LineWidth',1);
                hold on
%                 ax.YLim=[800,2400];
                ax.XLim=[-0.2,0.5];
                ax.TickDir='out';
                ax.Box='off';
                ax.TickLength=[0.02,0.02];
                ax.FontSize=6;
%                 ax.YTickLabel={};
            end
        end
    end   
end
hold on
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
filename=['_ori',num2str(id_ori),];
exportgraphics(gcf,[mysaveplotPath,'/raster/',filename,'_vector.pdf'],'ContentType','vector') 
close all           
%%
% for i= 1:Ncell_included_all %80% 91%[98,99,102,115,131,140]% 90% 86% 60 %, example[98,96,86,60]
for i=[122,123,126] % i= 1:Ncell_included_all %80% 91%[98,99,102,115,131,140]% 90% 86% 60 %, example[98,96,86,60]
    sptime=sp.st(sp.clu==cluster.depthsorted_id_includedall(i)); 
%     figure('Color',[1 1 1],'WindowState', 'maximized');    
    figure('Color',[1 1 1],'Position',[100 100 245 70],'Visible','on');   
    id_sz=1;
    for id_side=1:2
        for id_ori=2
            for id_lc=1:2
                %% PSTH of firing rate
                ax=subplot(3,4,[4+id_lc,8+id_lc]);                 
                p=plot(bincenter_psth-delay_sys,squeeze(psth_cond_included_sz_all([1,2]+(id_lc-1)*2,id_ori,:,i,id_sz)),'LineStyle','-','Linewidth',1);
                p(1).Color=colorlabel_side(1,:);
                p(2).Color=colorlabel_side(2,:);                               
%                 ax.XMinorTick='on';
                ax.TickDir='out';
                ax.TickLength=[0.04,0.04];
                ax.Box='off';
                ylim=get(gca,'YLim');
%                 ax.YLim=[0,ylim(2)];
                ax.YLim=[0,ceil(0.1+max(psth_cond_included_sz_all(:,id_ori,:,i,id_sz),[],'all')./5)*5];
                ax.XLim=[-0.1,0.5];
                ax.FontSize=6;
                hold on                             
                %% bar plot of firing rate
                ax2=subplot(3,4,[7,11]);                 
%                 b=bar(FR_sz_all([1,2]+(id_lc-1)*2,id_ori,i,id_sz));
                xx=[mean(FR_mat_sz_all([1,2]+(1-1)*2,id_ori,:,i,id_sz),3);mean(FR_mat_sz_all([1,2]+(2-1)*2,id_ori,:,i,id_sz),3)];
                xx_error=[std(FR_mat_sz_all([1,2]+(1-1)*2,id_ori,:,i,id_sz),0,3)./sqrt(N_repetition);std(FR_mat_sz_all([1,2]+(2-1)*2,id_ori,:,i,id_sz),0,3)./sqrt(N_repetition)];
                b=bar(xx);
                hold on
                errorbar(1:4,xx,xx_error,'Color',[0,0,0],'LineStyle','none','CapSize',1,'LineWidth',0.3)
                hold on
                b.EdgeColor="none";
                b.FaceColor='flat';
                b.CData=[colorlabel_side;colorlabel_side];
                if id_ori==1 && id_lc==1 && id_sz==1
                    ax2.YLabel.String='FR (Hz)';
                elseif (id_ori+id_lc)>2
%                     ax2.YTickLabel={};
                end   
                ax2.XTickLabel={};
                aatemp=max(FR_sz_all(:,:,i,id_sz),[],'all');
                if aatemp>0
                    ax2.YLim=[0,ceil(aatemp./5)*5];
                else
                    ax2.YLim=[0,1*5];
                end
                ax2.TickDir='out';
                xx=get(gca,'YLim');
                ax2.YTick=0:5:xx(2);
                ax2.TickLength=[0.04,0.04];
                ax2.Box='off';
                ax2.FontSize=6;
                %% raster plot
                event_on_time_currentcond_side1=event_on_time(trialmat(:,c_ori)==id_ori&trialmat(:,c_side)==id_side&trialmat(:,c_lc)==id_lc&trialmat(:,c_sz)==3-id_sz); 
                [~, ~, rasterx_side1, rasterY_side1, ~, ~] = psthAndBA(sptime, event_on_time_currentcond_side1, [0,1], psth_binSize);
                ax3=subplot(3,4,(id_lc-1)*2+id_side);                 
                plot(rasterx_side1-delay_sys,rasterY_side1,'Color',colorlabel_side(id_side,:),'LineWidth',0.5);
                hold on
                ax3.YLim=[1,N_repetition+1];
                ax3.XLim=[-0.1,0.5];
%                 ax3.XLim=[0,1];
                ax3.TickDir='out';
                ax3.Box='off';
                ax3.TickLength=[0.04,0.04];
                ax3.YAxis.Visible='off';
                ax3.XTickLabel={};
                ax3.FontSize=6;
            end
        end   
    end
    hold on
    set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
    filename=[num2str(i),'_clu',num2str(cluster.depthsorted_id_includedall(i))];
    exportgraphics(gcf,[mysaveplotPath,'/raster/',filename,'_vector_s.pdf'],'ContentType','vector') 
%     close all           
end 




%%
% for i= 1:Ncell_included_all %80% 91%[98,99,102,115,131,140]% 90% 86% 60 %, example[98,96,86,60]
%     sptime=sp.st(sp.clu==cluster.depthsorted_id_includedall(i)); 
% %     figure('Color',[1 1 1],'WindowState', 'maximized');     
%     figure('Color',[1 1 1],'Position',[100 100 1600 1200],'Visible','on');     
%     for id_sz=1:2
%         for id_ori=1:4
%             for id_lc=1:2
%                 %% PSTH of firing rate
%                 ax=subplot(6,8,(id_sz-1)*24+(id_ori-1)*2+id_lc);                 
%                 p=plot(bincenter_psth-delay_sys,squeeze(psth_cond_included_sz_all([1,2]+(id_lc-1)*2,id_ori,:,i,id_sz)),'LineStyle','-','Linewidth',1);
%                 p(1).Color=colorlabel_side(1,:);
%                 p(2).Color=colorlabel_side(2,:);                               
% %                 ax.XMinorTick='on';
%                 ax.TickDir='out';
%                 ax.TickLength=[0.04,0.04];
%                 ax.Box='off';
%                 ylim=get(gca,'YLim');
% %                 ax.YLim=[0,ylim(2)];
%                 ax.YLim=[0,ceil(0.1+max(psth_cond_included_sz_all(:,id_ori,:,i,id_sz),[],'all')./5)*5];
%                 ax.XLim=[-0.2,0.5];
%                 ax.FontSize=6;
%                 hold on
%                 if id_ori==1 && id_lc==1 && id_sz==1
%                     ax.YLabel.String='FR (Hz)';
%                     title(['Ori\_',num2str(square_ori(id_ori)),', LC',num2str(id_lc)]);
%                 elseif (id_ori+id_lc)>2
% %                     ax.YTickLabel={};
%                 end
%                 if id_sz==1 && (id_ori+id_lc)>2
%                     title([num2str(square_ori(id_ori)),', LC',num2str(id_lc)]);                    
%                 end
%                 %% bar plot of firing rate
%                 ax2=subplot(6,8,(id_sz-1)*24+(id_ori-1)*2+id_lc+8);  
% %                 b=bar(FR_sz_all([1,2]+(id_lc-1)*2,id_ori,i,id_sz));
%                 xx=mean(FR_mat_sz_all([1,2]+(id_lc-1)*2,id_ori,:,i,id_sz),3);
%                 xx_error=std(FR_mat_sz_all([1,2]+(id_lc-1)*2,id_ori,:,i,id_sz),0,3)./sqrt(N_repetition);
%                 b=bar(xx);
%                 hold on
%                 errorbar(1:2,xx,xx_error,'Color',[0,0,0],'LineStyle','none')
%                 hold on
%                 b.FaceColor='flat';
%                 b.CData=colorlabel_side;
%                 if id_ori==1 && id_lc==1 && id_sz==1
%                     ax2.YLabel.String='FR (Hz)';
%                 elseif (id_ori+id_lc)>2
% %                     ax2.YTickLabel={};
%                 end   
%                 ax2.XTickLabel={};
%                 aatemp=max(FR_sz_all(:,:,i,id_sz),[],'all');
%                 if aatemp>0
%                     ax2.YLim=[0,ceil(aatemp./5)*5];
%                 else
%                     ax2.YLim=[0,1*5];
%                 end
%                 ax2.TickDir='out';
%                 xx=get(gca,'YLim');
%                 ax2.YTick=0:5:xx(2);
%                 ax2.TickLength=[0.04,0.04];
%                 ax2.Box='off';
%                 ax2.FontSize=6;
%                 if p_anova_includedall(i,id_ori,1)<0.05/4 || p_anova_includedall(i,id_ori,3)<0.05/4
%                     color_p=[0 1 0];
%                 else 
%                     color_p=[0 0 0];
%                 end
%                 if id_sz==1 && (id_ori+id_lc)>=2
%                     title([num2str(p_anova_includedall(i,id_ori,1),3),',',num2str(p_anova_includedall(i,id_ori,3),3)],'Color',color_p);                     
%                 end
%                 %% raster plot
%                 event_on_time_currentcond_side1=event_on_time(trialmat(:,c_ori)==id_ori&trialmat(:,c_side)==1&trialmat(:,c_lc)==id_lc&trialmat(:,c_sz)==3-id_sz); 
%                 [~, ~, rasterx_side1, rasterY_side1, ~, ~] = psthAndBA(sptime, event_on_time_currentcond_side1, [0,1], psth_binSize);
%                 event_on_time_currentcond_side2=event_on_time(trialmat(:,c_ori)==id_ori&trialmat(:,c_side)==2&trialmat(:,c_lc)==id_lc&trialmat(:,c_sz)==3-id_sz); 
%                 [~, ~, rasterx_side2, rasterY_side2, ~, ~] = psthAndBA(sptime, event_on_time_currentcond_side2, [0,1], psth_binSize);
%                 ax3=subplot(6,8,(id_sz-1)*24+(id_ori-1)*2+id_lc+16);                 
%                 plot(rasterx_side2-delay_sys,rasterY_side2,'Color',colorlabel_side(2,:),'LineWidth',0.5);
%                 hold on
%                 plot(rasterx_side1-delay_sys,rasterY_side1+N_repetition,'Color',colorlabel_side(1,:),'LineWidth',0.5);
%                 hold on
%                 ax3.YLim=[0,N_repetition*2];
%                 ax3.XLim=[-0.2,0.5];
% %                 ax3.XLim=[0,1];
%                 ax3.TickDir='out';
%                 ax3.Box='off';
%                 ax3.TickLength=[0.04,0.04];
%                 if id_ori==1 && id_lc==1 && id_sz==1
%                     ax3.YLabel.String='Size 8';
%                 elseif id_ori==1 && id_lc==1 && id_sz==2
%                     ax3.YLabel.String='Size 4';
%                 else
%                     ax3.YTickLabel={};
%                 end
%                 ax3.FontSize=6;
%             end
%         end   
%     end
%     suplabel(['Cluster=',num2str(cluster.depthsorted_id_includedall(i)),', ',unitlabel{cluster.depthsorted_label_includedall(i)},', Depth=',num2str(round(cluster.depthsort_includedall(i))),'\mum',',',' \color[rgb]{1 0 0}Sd1, \color[rgb]{0 0 1}Sd2'],'t',[0.08 0.08 0.88 0.88]);
%     hold on
%     set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 30 20])
%     filename=[num2str(i),'_clu',num2str(cluster.depthsorted_id_includedall(i))];
% %     exportgraphics(gcf,[mysaveplotPath,'/raster/',filename,'.png'],'resolution','300') 
%     exportgraphics(gcf,[mysaveplotPath,'/raster/',filename,'_vector.pdf'],'ContentType','vector') 
%     close all           
% end 






%%
% 
% for id_ori=1:4   
%     temp_max=max(psth_cond_included_sz8,[],[1,2,3],'omitnan');
%     psth_cond_included_sz8_norm=psth_cond_included_sz8./repmat(temp_max,[4,4,1000,1]);
%     psth_sz8=nan(2,size(psth_cond_included,2));
%     for id_lc=1:2
%         for id_side=1:2
%             im1=nexttile(id_ori+(id_side-1)*4+(id_lc-1)*8);
%             temp1=squeeze(psth_cond_included_sz8_norm(id_side+(id_lc-1)*2,id_ori,:,:)); 
%             psth_sz8(id_side,:)=mean(temp1');
%             image(bincenter_psth(1:1000),1:Ncell_included_all,temp1','CDataMapping','scaled');
%             colormap(im1,cmap)    
%             ax=gca;
%             if depth_isdeep==0
%                 ax.YDir='Normal';
%             end
%             ax.CLim=[0,1];
%             if id_ori>1
%                 ax.YAxis.Visible='off';
%             end     
%             if id_ori==4 && id_side==2
%                 colorbar
%             end 
%             if id_side==1
%                 ax.XAxis.Visible='off';
%             end 
%             hold on
%             xLimits = get(gca,'XLim');
%             if flag_exis_layer==1
%                 for temp_i=1:length(cell_layer_idx_border_included_i{1,id_ori})
%                     plot(xLimits,[cell_layer_idx_border_included_i{1,id_ori}(temp_i)+0.5,cell_layer_idx_border_included_i{1,id_ori}(temp_i)+0.5],'k--','LineWidth',layerline_width)
%                     hold on
%                     plot([0.1,0.1],[1,Ncell_included_all],'w--','LineWidth',layerline_width)
%                 end 
%             end
% 
%         end        
%     end
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     exportgraphics(gcf,[mysaveplotPath,'/','psth.pdf'],'ContentType','vector');  
