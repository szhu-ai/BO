N_waveform_template=size(tempPeakWF,1);
N_waveform_timepoints=size(tempPeakWF,2);
time_wf=1:N_waveform_timepoints;
dx=0.04;
dy=2;
NS_BS_boundary = 0.2 * 30;
MS_BS_boundary = 0.3 * 30;
[~,temp_trough]=min(tempPeakWF,[],2);
[~,temp_peak]=max(tempPeakWF,[],2);

celltype=nan(1,N_waveform_template);
celltype(tempDur<=NS_BS_boundary)=1;
celltype(tempDur>NS_BS_boundary&tempDur<=MS_BS_boundary)=2;
celltype(tempDur>MS_BS_boundary)=3;
celltype(temp_trough>temp_peak)=0;
aa=temp_peak-temp_trough;
tempDur(celltype'==0)=aa(celltype'==0);

waveform_colorlabel=[55 126 184;77 175 74;154,50,188;228 26 28]/255;
depth_layer_upper=[depth_L23(2);depth_L4b(2);depth_L4c(2);depth_L56(2);depth_WM(2)];
sz_font=8;
save(fullfile(myresultPath,recordingDate,recordingSession,'cell_type_index_shude.mat'),'celltype');
%% plot waveform
figure('Color',[1 1 1],'WindowState', 'normal');
for id_celltype=[0,3,2,1]
    plot(repmat(templateXpos(celltype==id_celltype)',N_waveform_timepoints,1)+dx*repmat(time_wf',1,length(find(celltype==id_celltype))),repmat(templateYpos(celltype==id_celltype)',N_waveform_timepoints,1)+dy*tempPeakWF(celltype==id_celltype,:)','Color',waveform_colorlabel(id_celltype+1,:),'LineWidth',2.5);
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
print('-vector','-dpdf', [mysaveplotPath,'\',recordingDate,recordingSession,'_waveform_',testtype,'.pdf'], '-r0');  


%% plot waveform hot map

depthsort_temp=cluster_depth'+rand(1,N_cluster);
depth_int=interp1(depthsort_temp,depthsort_temp,1:2500,'nearest');
depth_gap=find(diff(depthsort_temp)>=1000);
if ~isempty(depth_gap)
    for jjj=1:length(depth_gap)
        temp0=depthsort_temp(depth_gap(jjj));
        temp1=depthsort_temp(depth_gap(jjj)+1);
        depth_int(round(temp0+5):round(temp1-5))=nan;
    end
end
Nrows=ceil(max(depthsort_temp)./500)*500;
wf_cmap=zeros(Nrows,N_waveform_timepoints);
cluster_waveform_norm=cluster_peakWF./max(cluster_peakWF,[],2);
for jjj=1:N_cluster
    wf_cmap(depth_int==depthsort_temp(jjj),:)=repmat(cluster_waveform_norm(jjj,:),length(find(depth_int==depthsort_temp(jjj))),1);
end
if depth_isdeep~=1
    depth_maxselection=max(depth_L23(2));
end
wf_cmap_selected=wf_cmap(1:depth_maxselection,:);
% figure('Position',[400,400,800,800],'Units', 'points');
% ax_pos=[0.2, 0.2, 0.6, 0.6];
% if sessionidx==10 || sessionidx==13
%     ax_pos=[0.2, 0.2, 0.16, 0.2];
%     colorbar_pos=[0.368,0.2,0.01,0.2];   
% else
%     ax_pos=[0.2, 0.2, 0.16, 0.3];
%     colorbar_pos=[0.368,0.2,0.01,0.3];
% end

% ax=axes('Position', ax_pos,'Units', 'points');
figure;
im1=image(time_wf/30,(1:depth_maxselection),wf_cmap_selected,'CDataMapping','scaled');
colormap(jet)
ax=gca;
ax.CLim=[-2,2];
% colorbar(ax,'Ticks',[],'Box','off','Color','none','Position',colorbar_pos,'Units','points')

hold on
xLimits = get(gca,'XLim');

for id_layer=1:5
    plot(ax,xLimits,[depth_layer_upper(id_layer),depth_layer_upper(id_layer)],'Color',[0 0 0 0.5],'LineWidth',layerline_width,'LineStyle','--')
    hold on
end

% plot(ax,[0,0],[-1300,1500],'Color',[0 0 0],'LineWidth',1);
% hold on
ax=gca;
ax.YTick=[0:100:2500];
ax.YLim=[0,2500];
ax.YDir='Normal';
ax.TickDir = 'out';
ax.Box='off';
ax.FontWeight='Bold';
ax.FontSize=sz_font;    
ax.YAxis.LineWidth = 1;
ax.XAxis.LineWidth = 1;

%% plot distribution of waveform duration, bar graph for all and all layers
% clusterDur_t=tempDur/30;
% % clusterDur_selected=clusterDur_t(cell_layer_index>=0);
% clusterDur_selected=clusterDur_t;
% 
% edges1=-0.6:1/30:0.6;
% edges1=edges1+0.001;
% edges2=0:1/30:0.6;
% edges2=edges2+0.001;
% 
% figure('Position',[400,400,500,500],'Units', 'points');
% ax_pos=[0.1, 0.1, 0.8, 0.8];
% ax=axes('Position', ax_pos,'Units', 'points');
% histogram(ax,clusterDur_selected(clusterDur_selected>0&clusterDur_selected<=0.201),edges1,'FaceColor',waveform_colorlabel(2,:),'EdgeColor','none','FaceAlpha',1)
% hold on
% histogram(ax,clusterDur_selected(clusterDur_selected>0.201&clusterDur_selected<=0.301),edges1,'FaceColor',waveform_colorlabel(3,:),'EdgeColor','none','FaceAlpha',1)
% hold on
% histogram(ax,clusterDur_selected(clusterDur_selected>0.301),edges1,'FaceColor',waveform_colorlabel(4,:),'EdgeColor','none','FaceAlpha',1)
% hold on
% histogram(ax,clusterDur_selected(clusterDur_selected<=0),edges1,'FaceColor',waveform_colorlabel(1,:),'EdgeColor','none','FaceAlpha',1)
% ax=gca;
% ax.XLim=[-0.6,0.6];
% ax.TickDir = 'out';
% ax.Box='off';    
% ax.FontWeight='Bold';
% ax.FontSize=sz_font;    
% ax.YAxis.LineWidth = 1;
% ax.XAxis.LineWidth = 1;  
% % filename=['WFdur_dis_',num2str(N_celltype)];
% % print('-vector','-dpdf', [mysaveplotPath,'\',filename,'.pdf']); 