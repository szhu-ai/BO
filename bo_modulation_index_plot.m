%% histogram of BO modulation index
close all
ccmap=brewermap(11,'PiYG');
colorlabel_sig=[67/255,147/255,195/255;ccmap(end-1,:);ccmap(2,:);ccmap(2,:)];
bo_lable={'Side1','Side2','Sideall'};
ori_lable={'ori0','ori45','ori90','ori135'};
lc_lable={'LC1','LC2','LCall'};
x_edge=-1.05:0.1:1.05;
mysaveplotPath_mod=fullfile(mysaveplotPath,'Modulation index');
%% plot side modulation index histogram
for id_lc=[1,2,3]
    figure('Color',[1 1 1],'Position',[100 100 1000 130]);   
    sgtitle(['Side Modulation index, ',lc_lable{id_lc}])
    for id_ori=1:4
        for id_layer=0
            % subplot(5,4,(4-id_layer)*4+id_ori)
            subplot(1,4,id_ori)
            if id_layer==0
                BO_selected=BO_index_normbymax(Includedidx_all,id_ori,id_lc);
                histogram(BO_selected,x_edge,'FaceColor',[0.5 0.5 0.5],'EdgeColor','none')
            else
                BO_selected=BO_index_normbymax(Includedidx_all&cell_layer_idx'==id_layer,id_ori,id_lc);
                histogram(BO_selected,x_edge,'FaceColor',colorlabel_layer(id_layer,:))
            end
            hold on
            yy=get(gca,'YLim');
            plot([nanmedian(BO_selected),nanmedian(BO_selected)],[0,yy(2)+2],'r--','LineWidth',0.5)   
            try
                ci=bootci(1000,@median,BO_selected);
            catch
                ci=[nan,nan];
            end
            hold on
            plot([0,0],[0,yy(2)+2],'k--','LineWidth',0.5)
            hold on
            ax=gca;
            ax.XLim=[-1.2,1.2];
            ax.TickDir = 'out';
            ax.Box='off';
            ax.TickLength=[0.02,0.02];
            ax.FontSize=12;
            try
                p=signtest(BO_selected);
            catch
                p=nan;
            end
            if id_ori>=1
                str_temp=[', Md=',num2str(nanmedian(BO_selected),3)];
            else
                str_temp='';
            end
            if p<0.05
                title({['p=',num2str(p,2),str_temp],[num2str(ci(1),2),'\_',num2str(ci(2),2)]})
            else
                title({['p=',num2str(p,2),str_temp],[num2str(ci(1),2),'\_',num2str(ci(2),2)]})                
            end
        end
    end
    exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_histo_',lc_lable{id_lc},'.pdf'],'ContentType','vector');  
end
%% correlation of side modulation 
% figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
% sgtitle(['Side Modulation index, correlation'])
% for id_ori=1:4
%     subplot(1,4,id_ori)    
%     for id_layer=1:4
%         BO1_selected=BO_index(Includedidx_all&cell_layer_idx'==id_layer,id_ori,1);
%         BO2_selected=BO_index(Includedidx_all&cell_layer_idx'==id_layer,id_ori,2);
%         scatter(BO1_selected,BO2_selected,20,colorlabel_layer(id_layer,:),'filled');
%         hold on
%         xlim([-1,1])
%         ylim([-1,1])
%         axis square
%         if id_ori==1
%             xlabel('Side modulation\_LC1')
%             ylabel('Side modulation\_LC2')
%         end
%     end
% end
% exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_corr.pdf'],'ContentType','vector'); 

%% depth plot of BO modulation index, colorcoded by anova_side only 
% close all
% figure('Color',[1 1 1],'Position',[100 100 550 1000]);   
% sgtitle(['Side Modulation index, ','N=',num2str(Ncell_included_all)])
% N_side_sig=nan(4,4,2);
% N_side_sig_anylayer=nan(2,4);
% p_anova_side_includedall_anyori_anylc=min(p_anova_side_includedall,[],[2,3]);
% N_side_sig_anyori_anylayer_anylc=length(find(p_anova_side_includedall_anyori_anylc<0.05))./length(find(Includedidx_all));
% 
% for id_ori=1:4
%     N_side_sig_anylayer(1,id_ori)=length(find(p_anova_side_includedall(:,id_ori,1)<0.05))./length(find(Includedidx_all));
%     N_side_sig_anylayer(2,id_ori)=length(find(p_anova_side_includedall(:,id_ori,2)<0.05))./length(find(Includedidx_all));
% end
% for id_ori=1:4
%         for id_layer=1:4
%             N_side_sig(id_layer,id_ori,1)=length(find(p_anova_side_includedall(:,id_ori,1)<0.05/4&cell_layer_idx_includedall'==id_layer))./length(find(Includedidx_all&cell_layer_idx'==id_layer));
%             N_side_sig(id_layer,id_ori,2)=length(find(p_anova_side_includedall(:,id_ori,2)<0.05/4&cell_layer_idx_includedall'==id_layer))./length(find(Includedidx_all&cell_layer_idx'==id_layer));
%         end
%     for id_lc=1:2        
%         subplot(3,4,(id_lc-1)*4+id_ori)
%         if depth_isdeep==0
%             scatter(BO_index_included(p_anova_side_includedall(:,id_ori,id_lc)>=0.05,id_ori,id_lc),cluster.depthsort_includedall(p_anova_side_includedall(:,id_ori,id_lc)>=0.05)-depth_relativezero,20,[0.8 0.8 0.8],'filled');
%             hold on
%             scatter(BO_index_included(p_anova_side_includedall(:,id_ori,id_lc)<0.05,id_ori,id_lc),cluster.depthsort_includedall(p_anova_side_includedall(:,id_ori,id_lc)<0.05)-depth_relativezero,20,[0 0 0],'filled');
%             hold on
%         else
%             scatter(BO_index_included(p_anova_side_includedall(:,id_ori,id_lc)>=0.05,id_ori,id_lc),depth_relativezero-cluster.depthsort_includedall(p_anova_side_includedall(:,id_ori,id_lc)>=0.05),20,[0.8 0.8 0.8],'filled');
%             hold on
%             scatter(BO_index_included(p_anova_side_includedall(:,id_ori,id_lc)<0.05,id_ori,id_lc),depth_relativezero-cluster.depthsort_includedall(p_anova_side_includedall(:,id_ori,id_lc)<0.05),20,[0 0 0],'filled');
%             hold on
%         end
%         plot([0,0],[-1500,2500],'k--')
%         xlim([-1.2 1.2]);
%         if depth_isdeep==0
%             ylim([-800,1400]);
%         else
%             ylim([-500,900]);
%         end
%         xLimits = get(gca,'XLim');
%         if flag_exis_layer==1
%             for id_layer=1:4
%                 if depth_isdeep==0
%                     plot([-1.2,-1.2],[depth_mat(id_layer,1)-depth_relativezero,depth_mat(id_layer,2)-depth_relativezero],'Color',colorlabel_layer(id_layer,:),'LineWidth',layerline_width)
%                 else
%                     plot([-1.2,-1.2],[depth_relativezero-depth_mat(id_layer,2),depth_relativezero-depth_mat(id_layer,1)],'Color',colorlabel_layer(id_layer,:),'LineWidth',layerline_width)
%                 end
%                 hold on
%             end 
%         end
%         if id_ori==1
%             ylabel(lc_lable{id_lc});
%         end
%         ax=gca;
%         ax.TickDir = 'out';   
%         ax.TickLength=[0.02,0.02];
%         ax.YTick=-1000:500:2000;
%         ax.FontSize=7;
%         ax.Box='off';
%         hold on
%     end
% end
% exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_depth_anova_side','.pdf'],'ContentType','vector');  



%% depth plot of BO modulation index 
p_sig=0.05;
figure('Color',[1 1 1],'Position',[100 100 1000 360]);   
sgtitle(['Side Modulation index, ','N=',num2str(Ncell_included_all)])
N_BO_sig_layer=nan(4,4,4);
N_BO_sig=nan(4,4);
Ncell_layer=nan(1,4);
for id_ori=1:4
       id_nonsig=p_anova_includedall(:,id_ori,1)>=p_sig&p_anova_includedall(:,id_ori,2)>=p_sig;
    id_sig_type1=p_anova_includedall(:,id_ori,1)< p_sig&p_anova_includedall(:,id_ori,3)>=p_sig;
    id_sig_type2=p_anova_includedall(:,id_ori,1)>=p_sig&p_anova_includedall(:,id_ori,2)< p_sig;
    id_sig_type3=p_anova_includedall(:,id_ori,1)< p_sig&p_anova_includedall(:,id_ori,3)< p_sig;
    N_BO_sig(id_ori,1)=length(find(id_sig_type1));
    N_BO_sig(id_ori,2)=length(find(id_sig_type2));
    N_BO_sig(id_ori,3)=length(find(id_sig_type3));
    N_BO_sig(id_ori,4)=length(find(id_nonsig));

    for id_layer=1:4
        Ncell_layer(id_layer)=length(find(Includedidx_all&cell_layer_idx'==id_layer));
        N_BO_sig_layer(id_layer,1,id_ori)=length(find(id_sig_type1&cell_layer_idx_includedall'==id_layer));
        N_BO_sig_layer(id_layer,2,id_ori)=length(find(id_sig_type2&cell_layer_idx_includedall'==id_layer));
        N_BO_sig_layer(id_layer,3,id_ori)=length(find(id_sig_type3&cell_layer_idx_includedall'==id_layer));
        N_BO_sig_layer(id_layer,4,id_ori)=length(find(id_nonsig&cell_layer_idx_includedall'==id_layer));
    end
    sz_scatter=40;
    for id_lc=3        
        % subplot(3,4,(id_lc-1)*4+id_ori)
        subplot(1,4,id_ori)
        if depth_isdeep==0
            scatter(BO_index_included_normbymax(id_nonsig,id_ori,id_lc),cluster.depthsort_includedall(id_nonsig)-depth_relativezero,sz_scatter,[0.8 0.8 0.8],'filled');
            hold on
            scatter(BO_index_included_normbymax(id_sig_type2,id_ori,id_lc),cluster.depthsort_includedall(id_sig_type2)-depth_relativezero,sz_scatter,[0.8 0.8 0.8],'filled');
            hold on
            scatter(BO_index_included_normbymax(id_sig_type3,id_ori,id_lc),cluster.depthsort_includedall(id_sig_type3)-depth_relativezero,sz_scatter,[0 0 0],'filled');
            hold on
            scatter(BO_index_included_normbymax(id_sig_type1,id_ori,id_lc),cluster.depthsort_includedall(id_sig_type1)-depth_relativezero,sz_scatter,[0 0 0],'filled');
            hold on            
        else
            scatter(BO_index_included_normbymax(id_nonsig,id_ori,id_lc),depth_relativezero-cluster.depthsort_includedall(id_nonsig),sz_scatter,[0.8 0.8 0.8],'filled');
            hold on
            scatter(BO_index_included_normbymax(id_sig_type2,id_ori,id_lc),depth_relativezero-cluster.depthsort_includedall(id_sig_type2),sz_scatter,[0.8 0.8 0.8],'filled');
            hold on
            scatter(BO_index_included_normbymax(id_sig_type3,id_ori,id_lc),depth_relativezero-cluster.depthsort_includedall(id_sig_type3),sz_scatter,[0 0 0],'filled');
            hold on
            scatter(BO_index_included_normbymax(id_sig_type1,id_ori,id_lc),depth_relativezero-cluster.depthsort_includedall(id_sig_type1),sz_scatter,[0 0 0],'filled');
            hold on 
        end
        plot([0,0],[-1500,2500],'k--','Linewidth',1)
        xlim([-1.2 1.2]);
        if depth_isdeep==0
            ylim([-800,1400]);
        else
            ylim([-500,900]);
        end
        xLimits = get(gca,'XLim');
        if flag_exis_layer==1
            for id_layer=1:4
                if depth_isdeep==0
                    plot([-1.2,-1.2],[depth_mat(id_layer,1)-depth_relativezero,depth_mat(id_layer,2)-depth_relativezero],'Color',colorlabel_layer(id_layer,:),'LineWidth',layerline_width)
                else
                    plot([-1.2,-1.2],[depth_relativezero-depth_mat(id_layer,2),depth_relativezero-depth_mat(id_layer,1)],'Color',colorlabel_layer(id_layer,:),'LineWidth',layerline_width)
                end
                hold on
            end 
        end
        if id_ori==1
            ylabel(lc_lable{id_lc});
        end
        ax=gca;
        ax.TickDir = 'out';   
        ax.TickLength=[0.02,0.02];
        ax.YTick=-1000:500:2000;
        ax.FontSize=12;
        ax.Box='off';
        hold on
    end
    % (length(find(id_sig_type1))+length(find(id_sig_type3)))./Ncell_included_all
end
exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_depth_p=0.05','.pdf'],'ContentType','vector');  

%% plot bo and ccg
% sz=100;
% p_sig=0.01;
% for id_ori=2
%            id_nonsig=p_anova_includedall(:,id_ori,1)>=p_sig&p_anova_includedall(:,id_ori,2)>=p_sig;
%         id_sig_type1=p_anova_includedall(:,id_ori,1)< p_sig&p_anova_includedall(:,id_ori,2)>=p_sig;
%         id_sig_type2=p_anova_includedall(:,id_ori,1)>=p_sig&p_anova_includedall(:,id_ori,2)< p_sig;
%         id_sig_type3=p_anova_includedall(:,id_ori,1)< p_sig&p_anova_includedall(:,id_ori,2)< p_sig;
%     for id_lc=3        
%         figure('Color',[1 1 1],'Position',[100 100 1000 1000]); 
%         for i=1:size(ccg_data.ccg_control,1)
%             id_pre=ccg_data.pre_id(i);
%             id_post=ccg_data.post_id(i);
%             id_pre_MI=BO_index_included(id_pre,id_ori,id_lc);
%             id_post_MI=BO_index_included(id_post,id_ori,id_lc);
%             if depth_isdeep==0
%                 id_pre_depth=cluster.depthsort_includedall(id_pre)-depth_relativezero;
%                 id_post_depth=cluster.depthsort_includedall(id_post)-depth_relativezero;
%             else
%                 id_pre_depth=depth_relativezero-cluster.depthsort_includedall(id_pre);
%                 id_post_depth=depth_relativezero-cluster.depthsort_includedall(id_post);
%             end            
%             if ccg_data.peak_lag(i)<-1 && (id_sig_type1(id_pre)|| id_sig_type3(id_pre)) && (id_sig_type1(id_post)|| id_sig_type3(id_post))  && ccg_data.pre_layerID(i)~=ccg_data.post_layerID(i) 
%                 plot([id_pre_MI,id_post_MI],[id_pre_depth,id_post_depth],'color',[0 0 1 0.3],'LineWidth',2);
%                 hold on
%                 scatter(id_pre_MI,id_pre_depth,sz,[0 0 0])
%                 hold on
%                 scatter(id_post_MI,id_post_depth,sz,[0 0 0])
%                 hold on
%             elseif ccg_data.peak_lag(i)>1 && (id_sig_type1(id_pre)|| id_sig_type3(id_pre)) && (id_sig_type1(id_post)|| id_sig_type3(id_post))  && ccg_data.pre_layerID(i)~=ccg_data.post_layerID(i) 
%                 plot([id_pre_MI,id_post_MI],[id_pre_depth,id_post_depth],'color',[1 0 0 0.3],'LineWidth',2);
%                 hold on
%                 scatter(id_pre_MI,id_pre_depth,sz,[0 0 0])
%                 hold on
%                 scatter(id_post_MI,id_post_depth,sz,[0 0 0])
%                 hold on
%             end             
%         end
%         plot([0,0],[-1500,2500],'k--')
%         xlim([-1.2 1.2]);
%         if depth_isdeep==0
%             ylim([-800,1400]);
%         else
%             ylim([-500,900]);
%         end
%         xLimits = get(gca,'XLim');
%         if flag_exis_layer==1
%             for id_layer=1:4
%                 if depth_isdeep==0
%                     plot([-1.2,-1.2],[depth_mat(id_layer,1)-depth_relativezero,depth_mat(id_layer,2)-depth_relativezero],'Color',colorlabel_layer(id_layer,:),'LineWidth',layerline_width)
%                 else
%                     plot([-1.2,-1.2],[depth_relativezero-depth_mat(id_layer,2),depth_relativezero-depth_mat(id_layer,1)],'Color',colorlabel_layer(id_layer,:),'LineWidth',layerline_width)
%                 end
%                 hold on
%             end 
%         end
%         if id_ori==1
%             ylabel(lc_lable{id_lc});
%         end
%         ax=gca;
%         ax.TickDir = 'out';   
%         ax.TickLength=[0.02,0.02];
%         ax.YTick=-1000:500:2000;
%         ax.FontSize=14;
%         ax.Box='off';
%         hold on
%         title(['lc=',num2str(id_lc),', ori=',num2str(id_ori)])
%     end
% end

 
%% depth plot of BO modulation index and separated by cell type
figure('Color',[1 1 1],'Position',[100,100,1400,900]);
sgtitle(['Side Modulation index, ',', N=',num2str(Ncell_included_all)])
for id_lc=1:3
    for id_celltype=1:3
        for id_ori=1:4
            subplot(3,12,(id_lc-1)*12+(id_celltype-1)*4+id_ori)
            id_sig_type1=cluster.depthsorted_celltype_includedall'==id_celltype&p_anova_includedall(:,id_ori,1)<0.05/4&p_anova_includedall(:,id_ori,2)>=0.05/4&p_anova_includedall(:,id_ori,3)>=0.05/4;
            id_sig_type2=cluster.depthsorted_celltype_includedall'==id_celltype&p_anova_includedall(:,id_ori,1)>=0.05/4&p_anova_includedall(:,id_ori,2)<0.05/4&p_anova_includedall(:,id_ori,3)>=0.05/4;
            id_sig_type3=cluster.depthsorted_celltype_includedall'==id_celltype&p_anova_includedall(:,id_ori,1)<0.05/4&p_anova_includedall(:,id_ori,2)<0.05/4;
            id_sig_type4=cluster.depthsorted_celltype_includedall'==id_celltype&p_anova_includedall(:,id_ori,3)<0.05/4&~id_sig_type3;
            id_nonsig=cluster.depthsorted_celltype_includedall'==id_celltype&p_anova_includedall(:,id_ori,1)>=0.05/4&p_anova_includedall(:,id_ori,2)>=0.05/4&p_anova_includedall(:,id_ori,3)>=0.05/4;
            scatter(BO_index_included(id_nonsig,id_ori,id_lc),find(id_nonsig),50,[0.8 0.8 0.8],'filled');
            hold on   
            scatter(BO_index_included(id_sig_type2,id_ori,id_lc),find(id_sig_type2),50,[0 0 0],'filled');
            hold on
            scatter(BO_index_included(id_sig_type4,id_ori,id_lc),find(id_sig_type4),50,[0 0.8 0],'filled');
            hold on
            scatter(BO_index_included(id_sig_type3,id_ori,id_lc),find(id_sig_type3),50,[0 0 0.8],'filled');
            hold on
            scatter(BO_index_included(id_sig_type1,id_ori,id_lc),find(id_sig_type1),50,[0.8 0 0],'filled');
            hold on
            plot([0,0],[0,Ncell_included_all],'k-')
            xlim([-1 1]);
            ylim([0,ceil(Ncell_included_all./5)*5])
            ax=gca;
            ax.TickDir = 'out';
            if id_ori>1
                ax.YAxis.Visible='off';
            end                
            hold on
            xLimits = get(gca,'XLim');
            if flag_exis_layer==1
                for temp_i=1:length(cell_layer_idx_border_included)
                    plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
                    hold on
                end 
            end
            if id_ori==1
                ylabel(lc_lable{id_lc});
            end
            if id_lc==1
                ax.Title.String=['\color[rgb]{0.8 0 0}',num2str(length(find(id_sig_type1))),',',...
                    '\color[rgb]{0 0 0}',num2str(length(find(id_sig_type2))),',',...
                    '\color[rgb]{0 0 0.8}',num2str(length(find(id_sig_type3))),',',...
                    '\color[rgb]{0 0.8 0}',num2str(length(find(id_sig_type4))),',',...
                    '\color[rgb]{0.8 0.8 0.8}',num2str(length(find(id_nonsig)))];
            end
        end
    end
end
exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_celltype','.pdf'],'ContentType','vector');  
%% correlation of side modulation and lc modulation
figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
sgtitle(['Side Modulation index, correlation'])
for id_ori=1:4
    subplot(1,4,id_ori)    
    for id_layer=1:4
        BO_selected=BO_index(Includedidx_all&cell_layer_idx'==id_layer,id_ori,3);
        LC_selected=LC_index(Includedidx_all&cell_layer_idx'==id_layer,id_ori,3);
        scatter(BO_selected,LC_selected,20,colorlabel_layer(id_layer,:),'filled');
        hold on
        xlim([-1,1])
        ylim([-1,1])
        axis square
        if id_ori==1
            xlabel('BO modulation')
            ylabel('LC modulation')
        end
    end
end
exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_LC_corr.pdf'],'ContentType','vector'); 
%%
figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
ori_pair=nchoosek(1:4,2);
for idx=1:6
    subplot(2,3,idx)    
    title(['BO modulation\_ori',num2str(ori_pair(idx,1)),' vs ',num2str(ori_pair(idx,2))])
    hold on
    for id_layer=1:4
        BO_selected1=BO_index_included(cell_layer_idx_includedall'==id_layer,ori_pair(idx,1),3);
        BO_selected2=BO_index_included(cell_layer_idx_includedall'==id_layer,ori_pair(idx,2),3);
        scatter(BO_selected1,BO_selected2,20,colorlabel_layer(id_layer,:),'filled');
        hold on
        xlim([-1,1])
        ylim([-1,1])
        axis square
        if id_ori==1
            xlabel('BO modulation')
            ylabel('BO modulation')
        end
    end
end
exportgraphics(gcf,[mysaveplotPath_mod,'/','BO_Ori_corr.pdf'],'ContentType','vector'); 

%% histogram of LC modulation index
x_edge=-1.05:0.1:1.05;
for id_bo=[1,2]
    figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
    sgtitle(['LC Modulation index, Side',num2str(id_bo)])
    for id_ori=1:4
        for id_layer=0:4
            subplot(5,4,(4-id_layer)*4+id_ori)
            if id_layer==0
                LC_selected=LC_index(Includedidx_all,id_ori,id_bo);
                histogram(LC_selected,x_edge,'FaceColor',[0.5 0.5 0.5])
            else
                LC_selected=LC_index(Includedidx_all&cell_layer_idx'==id_layer,id_ori,id_bo);
                histogram(LC_selected,x_edge,'FaceColor',colorlabel_layer(id_layer,:))
            end
            hold on
            yy=get(gca,'YLim');
            plot([median(LC_selected),median(LC_selected)],[0,yy(2)+2],'r-','LineWidth',1)        
            hold on
            plot([0,0],[0,yy(2)+2],'k-','LineWidth',1)
            hold on
            ax=gca;
            ax.TickDir = 'out';
            ax.Box='off';
            try
                p=signrank(LC_selected);
            catch
                p=nan;
            end
            if id_ori==1
                str_temp=[', N=',num2str(length(LC_selected))];
            else
                str_temp='';
            end
            if p<0.05
                title(['p=',num2str(p,2),str_temp],'Color','red')
            else
                title(['p=',num2str(p,2),str_temp])                
            end
        end
    end
    exportgraphics(gcf,[mysaveplotPath_mod,'/','LC index_histo_',bo_lable{id_bo},'.pdf'],'ContentType','vector');  
end
%%
%%%%%%%%%% depth plot of LC modulation index 
figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
sgtitle(['LC Modulation index, ','N=',num2str(Ncell_included_all)])

for id_ori=1:4
    id_sig_type1=p_anova_includedall(:,id_ori,1)<0.05/4&p_anova_includedall(:,id_ori,2)>=0.05/4&p_anova_includedall(:,id_ori,3)>=0.05/4;
    id_sig_type2=p_anova_includedall(:,id_ori,1)>=0.05/4&p_anova_includedall(:,id_ori,2)<0.05/4&p_anova_includedall(:,id_ori,3)>=0.05/4;
    id_sig_type3=p_anova_includedall(:,id_ori,1)<0.05/4&p_anova_includedall(:,id_ori,2)<0.05/4;
    id_sig_type4=p_anova_includedall(:,id_ori,3)<0.05/4&~id_sig_type3;
    id_nonsig=p_anova_includedall(:,id_ori,1)>=0.05/4&p_anova_includedall(:,id_ori,2)>=0.05/4&p_anova_includedall(:,id_ori,3)>=0.05/4;
    for id_bo=1:3
        subplot(3,4,(id_bo-1)*4+id_ori)
        plot([0,0],[0,Ncell_included_all],'k-')
        hold on
        scatter(LC_index_included(id_nonsig,id_ori,id_bo),find(id_nonsig),50,[0.8 0.8 0.8],'filled');
        hold on   
        scatter(LC_index_included(id_sig_type4,id_ori,id_bo),find(id_sig_type4),50,[0.8 0.8 0.8],'filled');
        hold on
        scatter(LC_index_included(id_sig_type3,id_ori,id_bo),find(id_sig_type3),50,[0.8 0.8 0.8],'filled');
        hold on
        scatter(LC_index_included(id_sig_type1,id_ori,id_bo),find(id_sig_type1),50,[0.8 0.8 0.8],'filled');
        hold on
        scatter(LC_index_included(id_sig_type2,id_ori,id_bo),find(id_sig_type2),50,[0 0 0],'filled');
        hold on
        xlim([-1 1]);
        ylim([0,ceil(Ncell_included_all./5)*5])
        ax=gca;
        ax.TickDir = 'out';
    %     if id_ori>1
    %         ax.YAxis.Visible='off';
    %     end                
        hold on
        xLimits = get(gca,'XLim');
        if flag_exis_layer==1
            for temp_i=1:length(cell_layer_idx_border_included)
                plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
                hold on
            end 
        end
        ax.Title.String=['\color[rgb]{0 0 0}',num2str(length(find(id_sig_type2))),', ',...
            '\color[rgb]{0.8 0.8 0.8}',num2str(length(find(id_nonsig))+length(find(id_sig_type4))+length(find(id_sig_type3))+length(find(id_sig_type1)))];
    end
end
exportgraphics(gcf,[mysaveplotPath_mod,'/','LC index_depth','.pdf'],'ContentType','vector');  