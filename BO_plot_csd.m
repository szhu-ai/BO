%%% load csd value and plot it in a hot map
clear all
clc
%%
rootpath='G:\npix\';
addpath(genpath(fullfile(rootpath,'code')));
testtype='BO';
myresultPath=fullfile(rootpath,'spike',testtype);

load(fullfile(rootpath,'result','bo_lut.mat'));
csddatapath='C:\Users\Shude\Dropbox\data_Moore_neuropixels\Neuropixels\ori\LFP_CSD\';
sessionidx=13;

csdname=ST.csd{sessionidx};
mysaveplotPath=[rootpath,'result\csd\'];
recordingDate=ST.recordingDate{sessionidx};
recordingSession=ST.recordingSession{sessionidx};
bankid=ST.bankid(sessionidx);

depth_maxselection=ST.ymax{sessionidx};
depth_isdeep=ST.isDeep(sessionidx);
depth_relativezero=ST.Zero{sessionidx};
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
    depth_L4b(1)=max(depth_L4b(1),0);
    if depth_L4b(1)>0
        depth_L23(1)=max(depth_L23(1),0);
    end
    depth_WM=[depth_L56(2),depth_maxselection];
end
chanelmapinfo=load('neuropixPhase3A_kilosortChanMap.mat');


%%
% Mcelltype=zeros(4,3);
% figure('Color',[1 1 1],'WindowState', 'maximized');
% for id_celltype=0:3
%     Mcelltype(1,id_celltype+1)=length(templateYpos(celltype'==id_celltype&templateYpos>=depth_L23(1)&templateYpos<depth_L23(2)));
%     Mcelltype(2,id_celltype+1)=length(templateYpos(celltype'==id_celltype&templateYpos>=depth_L4b(1)&templateYpos<depth_L4b(2)));
%     Mcelltype(3,id_celltype+1)=length(templateYpos(celltype'==id_celltype&templateYpos>=depth_L4c(1)&templateYpos<depth_L4c(2)));
%     Mcelltype(4,id_celltype+1)=length(templateYpos(celltype'==id_celltype&templateYpos>=depth_L56(1)&templateYpos<depth_L56(2)));
%     Mcelltype(5,id_celltype+1)=length(templateYpos(celltype'==id_celltype&templateYpos>=depth_WM(1)&templateYpos<depth_WM(2)));    
% end
% b=bar(Mcelltype);
% b(1).FaceColor=colorlabel(1,:);
% b(2).FaceColor=colorlabel(2,:);
% b(3).FaceColor=colorlabel(3,:);
% b(4).FaceColor=colorlabel(4,:);
% 
% %% plot CSD
% load([csddatapath,csdname,'.mat']);
% channeldepth=chanelmapinfo.ycoords;
% temp=reshape(channeldepth,4,[]);
% csd_depth=mean(temp,1);
% CSD_smooth=smoothdata(CSD,2,'lowess',25); %% 25*0.4ms=10ms
% csd_selected=CSD_smooth(csd_depth<=depth_maxselection,:);
% csd_selected_depth=csd_depth(csd_depth<=depth_maxselection);
% csd_clim=max(abs(csd_selected),[],'all');
% dy=400; %%1000;500;500;100;400
% figure('Color',[1 1 1],'WindowState', 'maximized');
% 
% im=image(time,csd_selected_depth,csd_selected,'CDataMapping','scaled');
% cmap=colormap(flipud(jet(256)));
% caxis([-csd_clim,csd_clim]);
% c=colorbar('Direction','reverse','Ticks',[],'Box','off','Color','none');
% hold on
% plot(time,csd_selected*dy+repmat(csd_selected_depth',1,size(time,2)),'Color',[0 0 0 0.5]);
% hold on
% layerline_width=2;
% for id_layer=1:5
%     plot([min(time),max(time)],[depth_layer_upper(id_layer),depth_layer_upper(id_layer)],'Color',[0 0 0 0.5],'LineWidth',layerline_width,'LineStyle','--')
%     hold on
% end
% ax=gca;
% ax.YLabel.String='Depth relative to Tip (\mum) ';
% ax.YDir='Normal';
% ax.TickDir = 'out';
% ax.PlotBoxAspectRatio = [1 1.5 1];
% ax.Box='off';
% ax.YLim=[0,depth_maxselection];
% print('-dpng', [mysaveplotPath,'\',csdname,'.png'], '-r0');
% print('-painters','-dpdf', [mysaveplotPath,'\',csdname,'.pdf'], '-r0');  