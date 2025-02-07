flag.stim_size='all';
flag.window_psth=[-0.3,1];
flag.start_time = 250; % ms, select spike train starting from 0 ms after stimulus onset
flag.end_time = 1000; % ms, select spike train until 1000 ms after stimulus onset
flag.start_time_index=flag.start_time-flag.window_psth(1)*1000; %data_subset = data.spiketrain(:, flag.start_time_index:flag.end_time_index, :);
flag.end_time_index=flag.end_time-flag.window_psth(1)*1000;
flag.jit_window = 25; % jitter spike timing within 25ms window, correlation slower than 25ms will be subtracted/corrected
flag.max_lag = 100; % maxium 100ms lag between two spike trains to computer CCG
flag.min_lag = -100; % better make it symmetrical for easy use purpose
flag.max_pt_lag = 10;   % select the center [-10,10]ms lag window for finding the location of peaks. 
flag.min_pt_lag = -10; 
flag.sig_max_lag = 10; % significance criteria, must have peak w/in 10 ms of 0

flag.sig_num_stds = 7; % must be 7 stds above noise mean
flag.sig_min_std = 0; % must have a non-zero noise std
% 1. data_for_ccg.spiketrain is a matrix with the size of N_trials*N_timebin*N_neurons
% 2. data_for_ccg.psth_window =[-0.3,1];
% 3. data_for_ccg.cluster_id is a vector of N_cluster*1, containing the id for each cluster
% 4. data_for_ccg.cluster_ypos is a vector of 1*N_cluster, containing the depth for each cluster

flag.Includedidx=Includedidx_all; % N_cluster*1
Ncell_included_ccg=length(find(flag.Includedidx));
flag.neuron_pair_all=nchoosek(1:Ncell_included_ccg,2);
flag.N_neuron_pair=size(flag.neuron_pair_all,1);
% get times of interest
if flag.stim_size==4
    id_cond_selected=find(rem(id_condition,2)==1);
    sz_label='_sz4';
elseif flag.stim_size==8
    id_cond_selected=find(rem(id_condition,2)==0);
    sz_label='_sz8';
else
    id_cond_selected=find(rem(id_condition,2)>=0);
    sz_label='_sz_all';   
end
%%
spiketrain_ccg=logical(data);
data_subset = spiketrain_ccg(id_cond_selected, flag.start_time_index:flag.end_time_index, flag.Includedidx);
data_subset_jitter=cell(Ncell_included_ccg,1);
data_subset_real=cell(Ncell_included_ccg,1);
% jitter spiketrains
% for j = 1:Ncell_included_ccg
%     data_subset_jitter{j} = BO_ccg_jitter(data_subset(:,:,j), flag.jit_window);
%     data_subset_real{j} = data_subset(:,:,j);
% end

for j = 1:Ncell_included_ccg
    data_subset_jitter{j,1} = zeros(size(data_subset(:,:,j)),'single');
    data_subset_real{j,1} = single(data_subset(:,:,j));
    for idx=1:N_condition
        data_subset_jitter{j,1}(id_cond_selected&id_condition==idx,:) = BO_ccg_jitter(data_subset(id_cond_selected&id_condition==idx,:,j), flag.jit_window);       
    end
end
%%
clearvars data_subset spiketrain_ccg
maxlag=flag.max_lag;
minlag=flag.min_lag;
n_pre=zeros(Ncell_included_ccg,maxlag-minlag+1,'single');
n_pre_jt=zeros(Ncell_included_ccg,maxlag-minlag+1,'single');

for id_neuron=1:Ncell_included_ccg
    for lag=-maxlag:1:-minlag
        if lag<=0
            prev = data_subset_real{id_neuron}(:,1:end+lag);
            prev_jt = data_subset_jitter{id_neuron}(:,1:end+lag);
        else
            prev = data_subset_real{id_neuron}(:,1+lag:end);
            prev_jt = data_subset_jitter{id_neuron}(:,1+lag:end);
        end
        n_pre(id_neuron,lag+maxlag+1)=sum(prev,'all');
        n_pre_jt(id_neuron,lag+maxlag+1)=sum(prev_jt,'all');         
    end
end
N_neuron_pair=flag.N_neuron_pair;
neuron_pair_all=flag.neuron_pair_all;
%%
ccgs=cell(N_neuron_pair,1);
parfor id_neuron_pair=1:N_neuron_pair
    disp("processing neuron pair: " + id_neuron_pair);
    ccgs{id_neuron_pair}=struct();
    pre_id=neuron_pair_all(id_neuron_pair,1);
    post_id=neuron_pair_all(id_neuron_pair,2);

    c_jk_norm = zeros(1,maxlag-minlag+1,'single');
    c_jk_unnorm= zeros(1,maxlag-minlag+1,'single');
    c_jk_norm_jt = zeros(1,maxlag-minlag+1,'single');
    c_jk_unnorm_jt= zeros(1,maxlag-minlag+1,'single');
    for lag=-maxlag:1:-minlag
        if n_pre(pre_id,lag+maxlag+1)*n_pre(post_id,-lag+maxlag+1)>0
            if lag<=0
                prev = data_subset_real{pre_id}(:,1:end+lag);
                postv = data_subset_real{post_id}(:,1-lag:end);
                prev_jt = data_subset_jitter{pre_id}(:,1:end+lag);
                postv_jt = data_subset_jitter{post_id}(:,1-lag:end);                
            else
                prev = data_subset_real{pre_id}(:,1+lag:end);
                postv = data_subset_real{post_id}(:,1:end-lag);
                prev_jt = data_subset_jitter{pre_id}(:,1+lag:end);
                postv_jt = data_subset_jitter{post_id}(:,1:end-lag);                
            end
            c_jk_unnorm(lag+maxlag+1) = sum(prev .* postv, 'all');
            c_jk_norm(lag+maxlag+1) = c_jk_unnorm(lag+maxlag+1)/sqrt(n_pre(pre_id,lag+maxlag+1)*n_pre(post_id,-lag+maxlag+1));
            c_jk_unnorm_jt(lag+maxlag+1) = sum(prev_jt .* postv_jt, 'all');
            c_jk_norm_jt(lag+maxlag+1) = c_jk_unnorm_jt(lag+maxlag+1)/sqrt(n_pre_jt(pre_id,lag+maxlag+1)*n_pre_jt(post_id,-lag+maxlag+1)); 
        end
    end
    ccgs{id_neuron_pair}.ccg_norm=c_jk_norm;
    ccgs{id_neuron_pair}.ccg_unnorm=c_jk_unnorm;
    ccgs{id_neuron_pair}.ccg_norm_jitter=c_jk_norm_jt;
    ccgs{id_neuron_pair}.ccg_unnorm_jitter=c_jk_unnorm_jt;
end

save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_ccg_data_',num2str(flag.start_time),'_',num2str(flag.end_time),sz_label,'_new.mat']),'ccgs','-v7.3');   
save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_ccg_flag_',num2str(flag.start_time),'_',num2str(flag.end_time),sz_label,'_new.mat']),'flag','-v7.3');   




%%
% spiketrain_ccg=logical(data);
% data_subset = spiketrain_ccg(id_cond_selected, flag.start_time_index:flag.end_time_index, flag.Includedidx);
% 
% % jitter spiketrains
% for j = 1:Ncell_included_ccg
%     data_subset_jitter{j} = BO_ccg_jitter(data_subset(:,:,j), flag.jit_window);
%     data_subset_real{j} = data_subset(:,:,j);
% end
% 
% clearvars data_subset spiketrain_ccg
% % change to for if do not want to use multiple cpus
% maxlag=flag.max_lag;
% minlag=flag.min_lag;
% 
% n_pre=zeros(Ncell_included_ccg,maxlag-minlag+1,'single');
% n_pre_jt=zeros(Ncell_included_ccg,maxlag-minlag+1,'single');
% 
% for id_neuron=1:Ncell_included_ccg
%     for lag=-maxlag:1:-minlag
%         if lag<=0
%             prev = data_subset_real{id_neuron}(:,1:end+lag);
%             prev_jt = data_subset_jitter{id_neuron}(:,1:end+lag);
%         else
%             prev = data_subset_real{id_neuron}(:,1+lag:end);
%             prev_jt = data_subset_jitter{id_neuron}(:,1+lag:end);
%         end
%         n_pre(id_neuron,lag+maxlag+1)=sum(prev,'all');
%         n_pre_jt(id_neuron,lag+maxlag+1)=sum(prev_jt,'all');         
%     end
% end
% N_neuron_pair=flag.N_neuron_pair;
% neuron_pair_all=flag.neuron_pair_all;
%%
% ccgs=cell(N_neuron_pair,1);
% for id_neuron_pair=560:600%1:N_neuron_pair
%     disp("processing neuron pair: " + id_neuron_pair);
%     ccgs{id_neuron_pair}=struct();
%     pre_id=neuron_pair_all(id_neuron_pair,1);
%     post_id=neuron_pair_all(id_neuron_pair,2);
% 
%     c_jk_norm = zeros(1,maxlag-minlag+1,'single');
%     c_jk_unnorm= zeros(1,maxlag-minlag+1,'single');
%     c_jk_norm_jt = zeros(1,maxlag-minlag+1,'single');
%     c_jk_unnorm_jt= zeros(1,maxlag-minlag+1,'single');
%     for lag=-maxlag:1:-minlag
%         if lag<=0
%             prev = data_subset_real{pre_id}(:,1:end+lag);
%             postv = data_subset_real{post_id}(:,1-lag:end);
%             prev_jt = data_subset_jitter{pre_id}(:,1:end+lag);
%             postv_jt = data_subset_jitter{post_id}(:,1-lag:end);                
%         else
%             prev = data_subset_real{pre_id}(:,1+lag:end);
%             postv = data_subset_real{post_id}(:,1:end-lag);
%             prev_jt = data_subset_jitter{pre_id}(:,1+lag:end);
%             postv_jt = data_subset_jitter{post_id}(:,1:end-lag);                
%         end
%         c_jk_unnorm(lag+maxlag+1) = sum(prev .* postv, 'all');
%         c_jk_norm(lag+maxlag+1) = c_jk_unnorm(lag+maxlag+1)/sqrt(n_pre(pre_id,lag+maxlag+1)*n_pre(post_id,-lag+maxlag+1));
%         c_jk_unnorm_jt(lag+maxlag+1) = sum(prev_jt .* postv_jt, 'all');
%         c_jk_norm_jt(lag+maxlag+1) = c_jk_unnorm_jt(lag+maxlag+1)/sqrt(n_pre_jt(pre_id,lag+maxlag+1)*n_pre_jt(post_id,-lag+maxlag+1));            
%     end
%     ccgs{id_neuron_pair}.ccg_norm=c_jk_norm;
%     ccgs{id_neuron_pair}.ccg_unnorm=c_jk_unnorm;
%     ccgs{id_neuron_pair}.ccg_norm_jitter=c_jk_norm_jt;
%     ccgs{id_neuron_pair}.ccg_unnorm_jitter=c_jk_unnorm_jt;
% end
% 
% figure;
% for id_neuron_pair=561:600
%     subplot(5,8,id_neuron_pair-560)
%     aa=ccgs{id_neuron_pair}.ccg_norm-ccgs{id_neuron_pair}.ccg_norm_jitter;
%     plot(-100:1:100,aa,'b')
% end


% ccgs=cell(N_neuron_pair,1);
% parfor id_neuron_pair=1:N_neuron_pair
%     disp("processing neuron pair: " + id_neuron_pair);
%     ccgs{id_neuron_pair}=struct();
%     pre_id=neuron_pair_all(id_neuron_pair,1);
%     post_id=neuron_pair_all(id_neuron_pair,2);
% 
%     c_jk_norm = zeros(1,maxlag-minlag+1,'single');
%     c_jk_unnorm= zeros(1,maxlag-minlag+1,'single');
%     c_jk_norm_jt = zeros(1,maxlag-minlag+1,'single');
%     c_jk_unnorm_jt= zeros(1,maxlag-minlag+1,'single');
%     for lag=-maxlag:1:-minlag
%         if lag<=0
%             prev = data_subset_real{pre_id}(:,1:end+lag);
%             postv = data_subset_real{post_id}(:,1-lag:end);
%             prev_jt = data_subset_jitter{pre_id}(:,1:end+lag);
%             postv_jt = data_subset_jitter{post_id}(:,1-lag:end);                
%         else
%             prev = data_subset_real{pre_id}(:,1+lag:end);
%             postv = data_subset_real{post_id}(:,1:end-lag);
%             prev_jt = data_subset_jitter{pre_id}(:,1+lag:end);
%             postv_jt = data_subset_jitter{post_id}(:,1:end-lag);                
%         end
%         c_jk_unnorm(lag+maxlag+1) = sum(prev .* postv, 'all');
%         c_jk_norm(lag+maxlag+1) = c_jk_unnorm(lag+maxlag+1)/sqrt(n_pre(pre_id,lag+maxlag+1)*n_pre(post_id,-lag+maxlag+1));
%         c_jk_unnorm_jt(lag+maxlag+1) = sum(prev_jt .* postv_jt, 'all');
%         c_jk_norm_jt(lag+maxlag+1) = c_jk_unnorm_jt(lag+maxlag+1)/sqrt(n_pre_jt(pre_id,lag+maxlag+1)*n_pre_jt(post_id,-lag+maxlag+1));            
%     end
%     ccgs{id_neuron_pair}.ccg_norm=c_jk_norm;
%     ccgs{id_neuron_pair}.ccg_unnorm=c_jk_unnorm;
%     ccgs{id_neuron_pair}.ccg_norm_jitter=c_jk_norm_jt;
%     ccgs{id_neuron_pair}.ccg_unnorm_jitter=c_jk_unnorm_jt;
% end
% %%
% save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_ccg_data_',num2str(flag.start_time),'_',num2str(flag.end_time),sz_label,'.mat']),'ccgs','-v7.3');   
% save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_ccg_flag_',num2str(flag.start_time),'_',num2str(flag.end_time),sz_label,'.mat']),'flag','-v7.3');   
