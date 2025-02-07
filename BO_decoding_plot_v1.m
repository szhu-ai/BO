%% set parameters for the timing counting window for different analysis
flag_trial_shuffle=2;
nameext_shuffle={'shuffle_','real_'};

flag_singleneuron=2;
nameext_neuron={'single_','population_'};

flag_singletime=1;
nameext_time={'onetime_','alltime_'};  

training_method=3; 
nameext_method={'classify_','svm_','ndt_'};  

if flag_singletime==1
    decoding_timebin=500;
    window_decoding=[0,0.5];
else
    decoding_timebin=50;
    window_decoding=[-0.1,0.5];
end
% 
load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_decording_error_',...
    nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},nameext_shuffle{2},num2str(decoding_timebin),'.mat']),'errs_all');   
if flag_singletime==1
    errs_shuffle=load(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_decording_error_',...
        nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},nameext_shuffle{1},num2str(decoding_timebin),'.mat']),'errs_all');   
end
window_decoding_idx=find(binborder_psth>=window_decoding(1),1,'first') :1: find(binborder_psth<=window_decoding(2),1,'last')-1;
bincenter_decoding=bincenter_psth(window_decoding_idx);
x_downsample=1:10:length(window_decoding_idx);
if flag_singletime==2
    bincenter_decoding_downsample=bincenter_decoding(x_downsample);
    N_timebin_decoding=length(x_downsample);
else
    N_timebin_decoding=1;
    bincenter_decoding_downsample=1;
    window_decoding=[0,0.5];
end

decode_pair_g=[1,2,3,4; 3,4,1,2; 1,3,2 4; 2,4,1,3];
N_decode_pair_g=size(decode_pair_g,1);

decode_pair_b=[1,2,3,4; 1,3,2 4];
N_decode_pair_b=size(decode_pair_b,1);

%%
N_neuron_decoding=size(errs_all,2);
n_ori=4;
cmap=brewermap(128,'*RdBu');
%%
ylable_string1={'Decode Side', 'Decode Side','Decode LC','Decode LC'};
ylable_string2={'train on LC1', 'train on LC2', 'train on Side1','train on Side2'};
if flag_singletime==1
    %%% plot average decoding performance using single time
    pf_generalize_real=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    pf_generalize_real_shuffle=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    pf_generalize_real_p=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);

    pf_generalize=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    pf_generalize_shuffle=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    pf_generalize_p=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    
    figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
    tiledlayout(4,4,'TileSpacing','tight');
    for id_decode_pair=1:4
        for id_ori=1:4
            for id_neurongroup=1:N_neuron_decoding
                x1=errs_all{1,id_neurongroup}.errs_generalize_real(id_decode_pair,:,:,id_ori);
                x2=errs_shuffle.errs_all{1,id_neurongroup}.errs_generalize_real(id_decode_pair,:,:,id_ori);
                pf_generalize_real(id_decode_pair,:,id_neurongroup,id_ori)=1-mean(x1,3);
                pf_generalize_real_shuffle(id_decode_pair,:,id_neurongroup,id_ori)=1-mean(x2,3);
%                 pf_generalize_real_p(id_decode_pair,:,id_neurongroup,id_ori)=ranksum(squeeze(x1),squeeze(x2));
                
                x1=errs_all{1,id_neurongroup}.errs_generalize(id_decode_pair,:,:,id_ori);
                x2=errs_shuffle.errs_all{1,id_neurongroup}.errs_generalize(id_decode_pair,:,:,id_ori);
                pf_generalize(id_decode_pair,:,id_neurongroup,id_ori)=1-mean(x1,3);
                pf_generalize_shuffle(id_decode_pair,:,id_neurongroup,id_ori)=1-mean(x2,3);
%                 pf_generalize_p(id_decode_pair,:,id_neurongroup,id_ori)=ranksum(squeeze(x1),squeeze(x2));
            end
            nexttile
            aa=squeeze(pf_generalize_real(id_decode_pair,:,:,id_ori));
            plot(aa*100,1:N_neuron_decoding,'r-','LineWidth',2);
            hold on
            ab=squeeze(pf_generalize(id_decode_pair,:,:,id_ori));
            plot(ab*100,1:N_neuron_decoding,'b-','LineWidth',2);
            hold on
%             scatter(0.05,find(pf_generalize_real_p(id_decode_pair,:,:,id_ori)>0.05/4),'green','filled')
            hold on
            plot([50,50],get(gca,'YLim'),'k--','LineWidth',1)
            xlim([0,100])
            ax=gca;
            ax.FontSize=16;
            ax.TickDir='out';
            ax.Box='off';
            hold on
            if id_decode_pair==1
                title(['Ori\_',num2str(id_ori)])
            elseif id_decode_pair==4
                xlabel('Performance (%)')
            end
            if id_ori>1
                ax.YAxis.Visible='off';
            else
                ylabel({ylable_string1{id_decode_pair},ylable_string2{id_decode_pair}})
            end      
            if flag_exis_layer==1
                xLimits = get(gca,'XLim');
                for temp_i=1:length(cell_layer_idx_border_included)
                    plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',1)
                    hold on
                end 
            end
        end
    end
    filename=['decoding_',nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},'generalized'];
    exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')
end
%%
ylable_string1={'Decode Side', 'Decode LC'};
ylable_string2={'train on all LC', 'train on all side'};
if flag_singletime==1
    pf_combined=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    pf_combined_shuffle=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    
    figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
    tiledlayout(2,4,'TileSpacing','tight');
    for id_decode_pair=1:2
        for id_ori=1:4
            for id_neurongroup=1:N_neuron_decoding
                N_temp=size(errs_all{1,id_neurongroup}.errs(id_decode_pair,:,:,id_ori),3);
                pf_combined(id_decode_pair,:,id_neurongroup,id_ori)=1-mean(errs_all{1,id_neurongroup}.errs(id_decode_pair,:,:,id_ori),3);
                pf_combined_shuffle(id_decode_pair,:,id_neurongroup,id_ori)=1-mean(errs_shuffle.errs_all{1,id_neurongroup}.errs(id_decode_pair,:,:,id_ori),3);
            end
            nexttile
            aa=squeeze(pf_combined(id_decode_pair,:,:,id_ori));
            plot(aa*100,1:N_neuron_decoding,'r-','LineWidth',2);
            hold on
            aa=squeeze(pf_combined_shuffle(id_decode_pair,:,:,id_ori));
            plot(aa*100,1:N_neuron_decoding,'k-','LineWidth',2);
            hold on
            plot([50,50],get(gca,'YLim'),'k--','LineWidth',1)
            ax=gca;
            ax.FontSize=16;
            ax.TickDir='out';
            ax.Box='off';
            hold on
            if id_ori>1
                ax.YAxis.Visible='off';
            else
                ylabel({ylable_string1{id_decode_pair},ylable_string2{id_decode_pair}})
            end    
            if id_decode_pair==1
                title(['Ori\_',num2str(id_ori)])
                xlim([20,80])
            elseif id_decode_pair==2
                xlabel('Performance (%)')
                xlim([0,100])
            end
            if flag_exis_layer==1
                xLimits = get(gca,'XLim');
                for temp_i=1:length(cell_layer_idx_border_included)
                    plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',1)
                    hold on
                end 
            end
        end
    end
    filename=['decoding_',nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime}];
    exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')
end
%% plot hotmap for population real
N_neuron_decoding=size(errs_all,2);
n_ori=4;
cmap=brewermap(128,'*RdBu');

ylable_string1={'Decode Side', 'Decode LC'};
ylable_string2={'train on all LC', 'train on all side'};
if flag_singletime==2
    pf_combined=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    pf_combined_std=nan(N_decode_pair_g,N_timebin_decoding,N_neuron_decoding,n_ori);
    
    figure('Color',[1 1 1],'Position',[100 100 800 1000]);   
    tiledlayout(2,4,'TileSpacing','tight');
    for id_decode_pair=1:2
        for id_ori=1:4
            for id_neurongroup=1:N_neuron_decoding
                N_temp=size(errs_all{1,id_neurongroup}.errs(id_decode_pair,:,:,id_ori),3);
                pf_combined(id_decode_pair,:,id_neurongroup,id_ori)=100-100*mean(errs_all{1,id_neurongroup}.errs(id_decode_pair,:,:,id_ori),3);
                pf_combined_std(id_decode_pair,:,id_neurongroup,id_ori)=100*std(1-errs_all{1,id_neurongroup}.errs(id_decode_pair,:,:,id_ori),[],3)./sqrt(N_temp);
    
            end
            nexttile
            aa=squeeze(pf_combined(id_decode_pair,:,:,id_ori));
            image(bincenter_decoding_downsample,1:N_neuron_decoding,aa','CDataMapping','scaled');
            colormap(cmap)
            xlim([-0.1,0.5])
            ax=gca;
            if depth_isdeep==0
                ax.YDir='Normal';
            end
            ax.CLim=[20,80];     
    %         ax.XTickLabel={'Ori0\newlineLC1','\newline2','Ori45\newline LC1','\newline2','Ori90\newline LC1','\newline2','Ori135\newline LC1','\newline2'};
            ax.FontSize=16;
            ax.TickDir='out';
            hold on
            if id_ori>1
                ax.YAxis.Visible='off';
            else
                ylabel({ylable_string1{id_decode_pair},ylable_string2{id_decode_pair}})
            end    
            if id_decode_pair==1
                title(['Ori\_',num2str(id_ori)])
            elseif id_decode_pair==2
                xlabel('Time (sec)')
            end
            if id_ori==4 && mod(id_decode_pair,2)==0
                colorbar
            end         
            if flag_exis_layer==1
                xLimits = get(gca,'XLim');
                for temp_i=1:length(cell_layer_idx_border_included)
                    plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
                    hold on
                end 
            end
        end
    end
    filename=['decoding_',nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},'psth_hot'];
    exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')
    %%%%%%%% plotted by condition
    figure('Color',[1 1 1],'Position',[100 100 1200 800]);   
    tiledlayout(2,4,'TileSpacing','tight');
    for id_decode_pair=1:2
        for id_ori=1:4
            nexttile
            for id_cell_layer=1:3
                pf_selected=squeeze(pf_combined(id_decode_pair,:,cluster.depthsorted_celllayer_includedall==id_cell_layer,id_ori));
                shadedErrorBar(bincenter_decoding_downsample,mean(pf_selected,2),std(pf_selected,[],2)./sqrt(size(pf_selected,2)),'lineProps',{'Color',colorlabel_layer(id_cell_layer,:),'LineWidth',2})
                hold on
            end
            plot(get(gca,'XLim'),[50,50],'k--','LineWidth',1)
            ax=gca;
            ax.FontSize=16;
            ax.TickDir='out';
            ax.XLim=[-0.1,0.5];
            ax.Box='off';
            ax.YLim=[30,80];
            hold on
            if id_ori>1
                ax.YAxis.Visible='off';
            else
                ylabel({ylable_string1{id_decode_pair},ylable_string2{id_decode_pair}})
            end    
            if id_decode_pair==1
                title(['Ori\_',num2str(id_ori)])
            elseif id_decode_pair==2
                xlabel('Time (sec)')
            end      
        end
    end
    filename=['decoding_',nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},'psth'];
    exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')
    %%%%%%% ploted by layer
    figure('Color',[1 1 1],'Position',[100 100 1200 800]);   
    for id_cell_layer=1:3
        for id_ori=1:4
            subplot(4,4,(4-id_cell_layer)*4+id_ori)
            for id_decode_pair=1:2
                pf_selected=squeeze(pf_combined(id_decode_pair,:,cluster.depthsorted_celllayer_includedall==id_cell_layer,id_ori));
                if id_decode_pair==1
                    shadedErrorBar(bincenter_decoding_downsample,mean(pf_selected,2),std(pf_selected,[],2)./sqrt(size(pf_selected,2)),'lineProps',{'Color',colorlabel_layer(id_cell_layer,:),'LineWidth',2})
                else
                    shadedErrorBar(bincenter_decoding_downsample,mean(pf_selected,2),std(pf_selected,[],2)./sqrt(size(pf_selected,2)),'lineProps',{'Color',[0,0,0],'LineWidth',2})
                end
                hold on
            end
            plot(get(gca,'XLim'),[50,50],'k--','LineWidth',1)
            ax=gca;
            ax.FontSize=16;
            ax.TickDir='out';
            ax.XLim=[-0.1,0.5];
            ax.Box='off';
            ax.YLim=[30,80];
            hold on
            if id_ori>1
                ax.YAxis.Visible='off';
            end
            if id_ori==1&&id_cell_layer==2
                ylabel('Decoding Performance (%)')
            end    
            if id_cell_layer==4
                title(['Ori\_',num2str(id_ori)])
            elseif id_cell_layer==1
                xlabel('Time (sec)')
            end        
        end
    end
    filename=['decoding_',nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},'psth_layer'];
    exportgraphics(gcf,[mysaveplotPath,'\',filename,'.pdf'],'ContentType','vector')
end
%%

%    
% 
%     if flag_exis_layer==3
%         bsl_pf_all=mean(pf_alltime(:,timebin_edges<=0.05,:,:),'all')+std(pf_alltime(:,timebin_edges<=0.05,:,:),0,'all');
%         sig_pf_all=nan(4,4,4); % (id_ori,id_decode_pair,id_layer)  
%         figure;    
%         tiledlayout(4,4,'TileSpacing','tight');
%         for id_decode_pair=1:4
%             for id_ori=1:4
%                 nexttile
%                 for id_layer=1:4
%                     pf_i=squeeze(pf_alltime(id_decode_pair,:,id_ori,cell_layer_idx==id_layer));
%                     h_pf=nan(size(pf_alltime,2),1);
%                     p_pf=h_pf;
%                     for id_timebin=1:size(pf_i,1)
%                         pf_temp=pf_i(id_timebin,:);
%     %                     [p_pf(id_timebin),h_pf(id_timebin)]=ranksum(pf_temp,0.5,'tail','right');
%                         [h_pf(id_timebin),p_pf(id_timebin)]=ttest(pf_temp,0.5,'tail','right');
%                     end
%                     xx=median(pf_i,2)>bsl_pf_all;
%                     xy=nan(length(xx)-1,1);
%                     for i=1:length(xx)-1
%                         xy(i)=(xx(i)&xx(i+1));
%                     end
%                     if mean(xy)>0
%                         sig_pf_all(id_ori,id_decode_pair,id_layer)=find(xy,1,'first');
%                     end
%                     plot(timebin_edges(1:end-1)+timebin_decoding/2,median(pf_i(timebin_edges(1:end-1)<=window_decoding(2),:),2),'Color',[colorlabel_layer(id_layer,:),0.9],'LineWidth',2)
%                     hold on
%                     ylim([0.4,1])
%                     xlim(window_decoding)
%                 end
%                 plot([0 0.5],[0.5 0.5],'--','LineWidth',2,'Color',[0 0 0])            
%                 ax=gca;
%                 ax.TickDir = 'out';
%                 ax.TickLength=[0.02,0.03];            
%                 ax.YTick=[0.4:0.1:1];
% %                 ax.XTick=[0:0.1:0.5];
% %                 ax.XTickLabel={'0','','','','','0.5'};            
%                 ax.YTickLabel={'','0.5','','','0.8','','1'};            
%                 ax.Box='off';
%                 if id_ori>1
%                     ax.YAxis.Visible='off';
%                 end
%                 if id_decode_pair<=5
%                     ax.XAxis.Visible='off';
%                 end
%                 if id_decode_pair<=2
%                     ax.YLim=[0.4,0.8];
%                 else
%                     ax.YLim=[0.4,1];
%                 end            
%             end
%         end
%     end



%     
% %     figure('Color',[1 1 1])%,'Position',[200,200,400,400],'Units', 'points');
%     pf_edges=0.05:0.1:1;
%     pf_reorder=pf([1,3,5,6,2,4],:,:);
%     pf_layer_median=nan(4,6,4,3);
%     for id_ori=1:4
%         for id_decode_pair=1:6
%             subplot(6,4,(id_decode_pair-1)*4+id_ori)
%             for id_layer=1:4
%                 pf_i=squeeze(pf_reorder(id_decode_pair,id_ori,cell_layer_idx==id_layer));
%                 [pf_N,edges]=histcounts(pf_i,pf_edges, 'Normalization', 'probability');
%                 pf_layer_median(id_ori,id_decode_pair,id_layer,:)=prctile(pf_i,[25,50,75]);
%                 hold on
%                 plot(repmat(id_layer,3,1),squeeze(pf_layer_median(id_ori,id_decode_pair,id_layer,:)),'-*','MarkerSize',5,'MarkerIndices',2,'Color',colorlabel_layer(id_layer,:),'LineWidth',2)
%                 hold on
%                 plot([0 4.5],[0.5 0.5],'--','LineWidth',2,'Color',[0 0 0])
%                 ylim([0.35 0.85])
%                 xlim([0.5 4.5])
%                 pf_cum=pf_N;
%                 for i=1:length(pf_N)
%                     pf_cum(i)=sum(pf_N(1:i));
%                 end
% %                 plot(pf_edges(1:end-1)+0.05,pf_N,'Color',colorlabel_layer(id_layer,:),'LineWidth',2)
% %                 hold on
%             end
% %             plot([0.5 0.5],[0 0.5],'--','LineWidth',2,'Color',[0 0 0])
% %             hold on
% %             plot([0 1],[0.5 0.5],'--','LineWidth',2,'Color',[0 0 0])            
%             ax=gca;
%             ax.TickDir = 'out';
%             ax.TickLength=[0.02,0.03];            
%             ax.YTick=[0.4:0.1:0.8];
%             ax.YTickLabel={'','0.5','','','0.8'};
%             
%             ax.Box='off';
%             if id_ori>1
%                 ax.YAxis.Visible='off';
%             end
%             if id_decode_pair<=6
%                 ax.XAxis.Visible='off';
%             end
%         end
% %         scatter(pf(id_ori,:),1:Ncell_included,100,'k','filled')
% %         xlim([0,1])
% %         for id_layer=1:4
% %             pd=fitdist(pf(id_ori,cell_layer_idx==id_layer)','gamma');
% %             x_values = 0:0.1:1;
% %             y = pdf(pd,x_values);
% %             plot(x_values,y,'LineWidth',2,'Color',colorlabel_layer(id_layer,:))
% %             hold on
% %         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    