clear all
close all
%%
rootpath='/Users/shushu/Dropbox/npix';
testtype='BO';
myresultPath=fullfile(rootpath,'spike',testtype);
load(fullfile(rootpath,'result','bo_lut.mat'));
mysaveplotPath=fullfile(rootpath,'result',testtype,'all');
%%
p_sig=0.01;
BO_MI_all=[];
BO_layer_index=[];
BO_depth=[];
BO_p_single=[];
BO_p_all=[];
BO_psth_all=[];
BO_sessionidx=[];
sessionidx=[5 6 7 10 11 13 15 16];
for idx=1:length(sessionidx)
    recordingDate=ST.recordingDate{sessionidx(idx)};
    recordingSession=ST.recordingSession{sessionidx(idx)};
    load(fullfile(myresultPath,'/all/',[recordingDate,'_',recordingSession,'_cell.mat']),'cell');
    BO_MI_all=[BO_MI_all,cell.BO_MI];
    BO_depth=[BO_depth,cell.depth];
    BO_layer_index=[BO_layer_index,cell.layer_index];
    BO_p_single=cat(3,BO_p_single,cell.p_single);
    BO_p_all=cat(3,BO_p_all,cell.p_all);
    BO_psth_all=cat(4,BO_psth_all,cell.psth); %% 4*4*500*N
    BO_sessionidx=[BO_sessionidx,repmat(sessionidx(idx),1,length(cell.layer_index))];
end
sysdelay=[70,70,70,70,70,70,1,1];
%%
BO_MI=reshape(BO_MI_all,4,2,[]);
BO_MI_double_sign=BO_MI(:,1,:).*BO_MI(:,2,:);
N_cell_layer=zeros(4,1);
N_sig_single=zeros(4,1);
N_sig_double=zeros(4,1);
N_sig_double_consistent=zeros(4,1);
N_sig_double_inconsistent=zeros(4,1);

N_sig_side=zeros(4,1);
N_sig_lc=zeros(4,1);
N_sig_interaction=zeros(4,1);
BO_MI_double_min_sign=nan(1,size(BO_MI,3));

p_single_min=transpose(squeeze(min(BO_p_single,[],[1 2]))); % 1*Ncell
p_double=max(BO_p_single,[],2);
[p_double_min_temp,idx_temp]=min(p_double);
for id_neuron=1:size(BO_MI,3)
    BO_MI_double_min_sign(id_neuron)=BO_MI_double_sign(idx_temp(id_neuron),:,id_neuron);
end
p_double_min=transpose(squeeze(p_double_min_temp));
p_side_min=transpose(squeeze(min(BO_p_all(:,1,:),[],1))); % 1*Ncell
p_lc_min=transpose(squeeze(min(BO_p_all(:,2,:),[],1))); % 1*Ncell
p_interaction_min=transpose(squeeze(min(BO_p_all(:,3,:),[],1))); % 1*Ncell
for id_layer=1:4
    N_cell_layer(id_layer)=length(find(BO_layer_index==id_layer));
    N_sig_single(id_layer)=length(find(p_single_min<p_sig&BO_layer_index==id_layer));
    N_sig_double(id_layer)=length(find(p_double_min<p_sig&BO_layer_index==id_layer));
    N_sig_double_consistent(id_layer)=length(find(p_double_min<p_sig&BO_layer_index==id_layer&BO_MI_double_min_sign>0));
    N_sig_double_inconsistent(id_layer)=length(find(p_double_min<p_sig&BO_layer_index==id_layer&BO_MI_double_min_sign<=0));    
    N_sig_side(id_layer)=length(find(p_side_min<p_sig&BO_layer_index==id_layer&p_lc_min>=p_sig));
    N_sig_lc(id_layer)=length(find(p_lc_min<p_sig&BO_layer_index==id_layer));
    N_sig_interaction(id_layer)=length(find(p_interaction_min<p_sig&BO_layer_index==id_layer));    
end
N_cell_total=sum(N_cell_layer);
prop_sig=N_sig_single./N_cell_layer;
prop_sig_consistent=N_sig_double_consistent./N_cell_layer;
prop_sig_inconsistent=N_sig_double_inconsistent./N_cell_layer;
N_sig_single-N_sig_double
%%
figure('Color',[1 1 1])
b=bar(flipud([N_cell_layer-N_sig_single,N_sig_single]));
b(1).FaceColor=[0 0 0];
ax=gca;
ax.Box='off';
ax.TickDir='out';
ax.YAxis.LineWidth = 1;
ax.XAxis.LineWidth = 1;
ax.TickLength=[0.02,0.03];
ax.FontSize=16;
exportgraphics(gcf,[mysaveplotPath,'\','counts.pdf'],'ContentType','vector')
%%
colorlabel_layer=brewermap(12,'Paired');
colorlabel_layer=colorlabel_layer([10,4,3,8,2],:);
BO_p_single_temp=nan(4,2,size(BO_MI,3));
BO_p_single_temp(:,1,:)=min(BO_p_single,[],2);
BO_p_single_temp(:,2,:)=min(BO_p_single,[],2);
aa=BO_MI.*(BO_p_single_temp<p_sig);
figure('Color',[1 1 1])
    scatter(reshape(squeeze(aa(:,1,BO_layer_index<=4)),1,[]),reshape(squeeze(aa(:,2,BO_layer_index<=4)),1,[]),50,'k','filled','MarkerFaceAlpha',0.6)
    axis square
ax=gca;
ax.XLim=[-1.2 1.2];
ax.XTick=[-1,0 1];
ax.YLim=[-1.2 1.2];
ax.YTick=[-1,0 1];
ax.TickDir='out';
ax.YAxis.LineWidth = 1;
ax.XAxis.LineWidth = 1;
ax.TickLength=[0.02,0.03];
ax.FontSize=16;

exportgraphics(gcf,[mysaveplotPath,'\','dependce on contrast.pdf'],'ContentType','vector')

% for idx_layer=1:4
%     scatter(reshape(squeeze(aa(:,1,BO_layer_index==idx_layer)),1,[]),reshape(squeeze(aa(:,2,BO_layer_index==idx_layer)),1,[]),50,'filled','MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.8)
%     hold on
% end
%%
prop_MI_consistent=nan(1,length(sessionidx));
BO_MI_sign=sign(BO_MI.*(BO_p_single<p_sig));
BO_MI_sign_temp=reshape(BO_MI_sign,8,[]);
x_sign=[];
x_dist=[];
x_session=[];
for idx=1:length(sessionidx)
    BO_MI_sign_temp_c=BO_MI_sign_temp(:,BO_sessionidx==sessionidx(idx)&BO_layer_index<=4&BO_layer_index>=1);
    BO_depth_c=BO_depth(:,BO_sessionidx==sessionidx(idx)&BO_layer_index<=4&BO_layer_index>=1);
    for idx_cond=1:8
        BO_MI_sign_temp_c_condi=BO_MI_sign_temp_c(idx_cond,:);
        if length(find(BO_MI_sign_temp_c_condi))>1
            c_vec=nchoosek(find(BO_MI_sign_temp_c_condi),2);
            xa=BO_MI_sign_temp_c_condi(c_vec);
            xb=BO_depth_c(c_vec);
            x_sign_temp=xa(:,1).*xa(:,2);
            x_dist_temp=xb(:,2)-xb(:,1);
            x_sign=[x_sign;x_sign_temp];
            x_dist=[x_dist;x_dist_temp];
            x_session=[x_session;repmat(sessionidx,size(c_vec,1),1)];
        end
    end
    a=sum(abs(BO_MI_sign(:,:,BO_sessionidx==sessionidx(idx)&BO_layer_index<=4&BO_layer_index>=1)),3);
    b=sum(BO_MI_sign(:,:,BO_sessionidx==sessionidx(idx)&BO_layer_index<=4&BO_layer_index>=1),3);
    b(a==1)=0;  %% exclude condiiton with only 1 sig neuron
    a(a==1)=0;
    c1=(a+b)/2;
    c2=(a-b)/2;
    d=max(c1,c2);
    prop_MI_consistent(idx)=sum(d)/sum(a);
end
figure('Color',[1 1 1])
tiledlayout(1,2,'TileSpacing','tight')
nexttile
scatter(1:length(sessionidx),prop_MI_consistent,100,'k','filled')
ax=gca;
ax.XLim=[0,length(sessionidx)+1];
ax.XTick=[1 3 5 7]+1;
ax.YLim=[0.4 1];
ax.TickDir='out';
ax.YAxis.LineWidth = 1;
ax.XAxis.LineWidth = 1;
ax.TickLength=[0.02,0.03];
ax.FontSize=16;
nexttile
edges=0:100:2000;
histogram(x_dist(x_sign==1),edges,'LineStyle','none')
hold on
histogram(x_dist(x_sign==-1),edges)
hold on
ax=gca;
ax.TickDir='out';
ax.Box='off';
ax.YAxis.LineWidth = 1;
ax.XAxis.LineWidth = 1;
ax.TickLength=[0.02,0.03];
ax.FontSize=16;
exportgraphics(gcf,[mysaveplotPath,'\','consistence.pdf'],'ContentType','vector')

%%  psth 4 cond          4 ori         500 bins        1229 cells
BO_psth_sig_side_layeridx=[];
BO_psth_sig_side_sessionidx=[];
BO_psth_sig_side=[];
for idx=1:length(sessionidx)
    for id_layer=1:4
        BO_psth_all_c=BO_psth_all(:,:,sysdelay(idx):sysdelay(idx)+300,BO_sessionidx==sessionidx(idx)&BO_layer_index==id_layer);
        BO_MI_sign_c=BO_MI_sign(:,:,BO_sessionidx==sessionidx(idx)&BO_layer_index==id_layer);
        if size(BO_MI_sign_c,3)>1
            for id_cell=1:size(BO_MI_sign_c,3)
                for id_ori=1:4
                    for id_lc=1:2
                        if BO_MI_sign_c(id_ori,id_lc,id_cell)==1
                          BO_psth_sig_side_temp= BO_psth_all_c([id_lc,id_lc+2],id_ori,:,id_cell);
                          BO_psth_sig_side=cat(4,BO_psth_sig_side,BO_psth_sig_side_temp);
                          BO_psth_sig_side_layeridx=[BO_psth_sig_side_layeridx,id_layer];
                          BO_psth_sig_side_sessionidx=[BO_psth_sig_side_sessionidx,sessionidx(idx)];                      
                        elseif BO_MI_sign_c(id_ori,id_lc,id_cell)==-1
                          BO_psth_sig_side_temp= BO_psth_all_c([id_lc+2,id_lc],id_ori,:,id_cell);
                          BO_psth_sig_side=cat(4,BO_psth_sig_side,BO_psth_sig_side_temp); 
                          BO_psth_sig_side_layeridx=[BO_psth_sig_side_layeridx,id_layer];
                          BO_psth_sig_side_sessionidx=[BO_psth_sig_side_sessionidx,sessionidx(idx)];                            
                        end
                    end
                end
            end
        end
    end
end
ccmap=brewermap(11,'PiYG');
colorlabel_sig=[67/255,147/255,195/255;ccmap(end-1,:);ccmap(2,:);ccmap(2,:)];
figure('Color',[1 1 1]);
tiledlayout(4,1,'TileSpacing','tight');
xtime=0.001:0.001:0.3;
for id_layer=4:-1:1
    nexttile
    yy=squeeze(BO_psth_sig_side(1,1,1:300,BO_psth_sig_side_layeridx==id_layer));
    shadedErrorBar(xtime,mean(yy,2),std(yy,0,2)./sqrt(size(yy,2)),'lineProps',{'Color',colorlabel_sig(3,:),'Linewidth',2});
    hold on
    yy=squeeze(BO_psth_sig_side(2,1,1:300,BO_psth_sig_side_layeridx==id_layer));
    shadedErrorBar(xtime,mean(yy,2),std(yy,0,2)./sqrt(size(yy,2)),'lineProps',{'Color',colorlabel_sig(2,:),'Linewidth',2});
    ax=gca;
    if id_layer>1
        ax.XAxis.Visible='off';
    end  
    ax.TickDir='out';
    ax.YLim=[0,0.6];
    ax.YAxis.LineWidth = 1;
    ax.XAxis.LineWidth = 1;
    ax.TickLength=[0.02,0.03];
    ax.FontSize=16;
end
exportgraphics(gcf,[mysaveplotPath,'\','psth.pdf'],'ContentType','vector')

% for idx=1:length(sessionidx)
%     figure('Color',[1 1 1]);
%     tiledlayout(4,1,'TileSpacing','tight');
%     for id_layer=4:-1:1
%         nexttile
%         plot(1:300,squeeze(mean(BO_psth_sig_side(1,1,sysdelay(idx):sysdelay(idx)+300,BO_psth_sig_side_layeridx==id_layer&BO_psth_sig_side_sessionidx==sessionidx(idx)),4)),'Color',colorlabel_sig(3,:),'LineWidth',2)
%         hold on
%         plot(1:300,squeeze(mean(BO_psth_sig_side(2,1,sysdelay(idx):sysdelay(idx)+300,BO_psth_sig_side_layeridx==id_layer&BO_psth_sig_side_sessionidx==sessionidx(idx)),4)),'Color',colorlabel_sig(2,:),'LineWidth',2)    
%     end
% end