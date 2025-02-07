close all
%% set parameters for the timing counting window for different analysis
flag_trial_shuffle=1;
nameext_shuffle={'shuffle_','real_'};

flag_singletime=1;
nameext_time={'onetime_','alltime_'};  

flag_singleneuron=2;
nameext_neuron={'single_','population_'};

training_method=3; 
nameext_method={'classify_','svm_','ndt_'};  

decoding_timebin=50;
if flag_singletime==1
    decoding_timebin=500;
    window_decoding=[0,0.5];
else
    window_decoding=[-0.1,0.5];
end
%%
toolbox_dir = '/Users/shushu/Dropbox/npix/code_new/decoding/ndt.1.0.4/';
addpath(toolbox_dir)
%% compute or load the firing rate for each condition
[~,sortidx]=sort(id_condition);
aa=reshape(data(sortidx,:,:),N_repetition,N_condition,N_Bins,N_cluster)./psth_binSize;
spikes_all_temp=permute(aa,[4,2,3,1]); %%(N_cluster,N_condition,N_timebin,N_repetition);
clear aa
spikes_all_temp=spikes_all_temp(:,2:2:32,:,:);  
temp_max_sz8=squeeze(max(psth_cond_all(2:2:32,window_spikecount_idx,:),[],[1,2],'omitnan'));

spikes_all_temp=spikes_all_temp./repmat(temp_max_sz8,1,16,N_Bins,N_repetition);
%%
% figure;
% for id_ori=1:4
%     subplot(1,4,id_ori)
%     aa=spikes_all_temp(Includedidx_all,(id_ori-1)*4+1:id_ori*4,:,:);
%     plot(bincenter_psth',movmean(squeeze(mean(aa(:,1,:,:)-aa(:,2,:,:)+aa(:,3,:,:)-aa(:,4,:,:),[1,2,4])),50),'r')
%     hold on
%     plot(bincenter_psth',movmean(squeeze(mean(aa(:,1,:,:)+aa(:,2,:,:)-aa(:,3,:,:)-aa(:,4,:,:),[1,2,4])),50),'k')    
%     hold on
%     ylim([-0.3,0.4])
% end
%%    
window_decoding_idx=find(binborder_psth>=window_decoding(1),1,'first') :1: find(binborder_psth<=window_decoding(2),1,'last')-1;
bincenter_decoding=bincenter_psth(window_decoding_idx);
x_downsample=1:10:length(window_decoding_idx);

spikes_decoding=spikes_all_temp(Includedidx_all,:,window_decoding_idx,:);
if flag_singletime==1
    spikes_downsample=mean(spikes_decoding,3);
    N_timebin_decoding=1;
    bincenter_decoding_downsample=1;
else  % smooth data within 100 ms and then downsample by 20 times
    spikes_sm=movmean(spikes_decoding,decoding_timebin,3);
    spikes_downsample=single(spikes_sm(:,:,x_downsample,:));
%     spikes_bin10_temp=reshape(spikes_sm,Ncell_included_all,N_condition,10,[],N_repetition);
%     spikes_downsample=squeeze(mean(spikes_bin10_temp,3));
    N_timebin_decoding=size(spikes_downsample,3);
    bincenter_decoding_downsample=bincenter_decoding(x_downsample);
end
clear spikes_sm data_bsl spikes_decoding spikes_all_temp FR_mat_reshape_sz8_i FR_mat_reshape_sz8 FR_mat_reshape
%     spikes_avg=zscore(spikes_avg,0,[2,3,4]);        
%     figure('Color',[1 1 1]);%     aa=mean(spikes_avg,4);%     imagesc(aa./max(aa,[],2)) %     colormap(jet)
%%
decode_pair_g=[1,2,3,4; 3,4,1,2; 1,3,2 4; 2,4,1,3];
N_decode_pair_g=size(decode_pair_g,1);

decode_pair_b=[1,2,3,4; 1,3,2 4];
N_decode_pair_b=size(decode_pair_b,1);

nTest=1;                             % leave one out %   nTest = ceil(nTrials * .2);  
nTrials=size(spikes_downsample,4);  %% 2 binary class compariasion
cvt1 = nchoosek(1:nTrials,nTest);  
ncvt1=size(cvt1,1);
cv_reps=ncvt1*ncvt1; 
cv_reps_b=2000;
%         warning('off','all')   
if flag_singleneuron==1
    n_per_group=1;
    N_neuron_group=Ncell_included_all;
else
    n_per_group=9;
    id_neuron_downsample=1:Ncell_included_all;
    N_neuron_group=length(id_neuron_downsample);
end
n_ori = 4;
errs_all=cell(1,N_neuron_group);
parfor id_neurongroup=1:N_neuron_group
        errs_all{id_neurongroup}=struct();

        errs=zeros(N_decode_pair_b,N_timebin_decoding,cv_reps_b,n_ori,'single')+0.5;
        errs_generalize_real=zeros(N_decode_pair_g,N_timebin_decoding,cv_reps_b,n_ori,'single')+0.5;
        errs_generalize=zeros(N_decode_pair_g,N_timebin_decoding,cv_reps_b,n_ori,'single')+0.5;

        fprintf('Group %d \n',id_neurongroup)
        for id_ori = 1:4
            if flag_singleneuron==1
                c_spikes_side=spikes_downsample(id_neurongroup,(id_ori-1)*4+[1,2,3,4],:,:);
            elseif flag_singleneuron==2
                c0=id_neuron_downsample(id_neurongroup)-4;
                c1=id_neuron_downsample(id_neurongroup)+4;
                c0=max(c0,1);
                c1=min(c1,N_neuron_group);
                c_spikes_side=spikes_downsample(c0:c1,(id_ori-1)*4+[1,2,3,4],:,:);
            end
            %%  for generalize
            for c_decode_pair=1:N_decode_pair_g %%% for generalize
                for id_timebine=1:N_timebin_decoding   
                    c_spikes1=squeeze(c_spikes_side(:,decode_pair_g(c_decode_pair,1),id_timebine,:));
                    c_spikes2=squeeze(c_spikes_side(:,decode_pair_g(c_decode_pair,2),id_timebine,:));
                    c_spikes3=squeeze(c_spikes_side(:,decode_pair_g(c_decode_pair,3),id_timebine,:));
                    c_spikes4=squeeze(c_spikes_side(:,decode_pair_g(c_decode_pair,4),id_timebine,:));                   
                    if ~isvector(c_spikes1)
                        c_spikes1 =c_spikes1';                        
                        c_spikes2 =c_spikes2';
                        c_spikes3 =c_spikes3';
                        c_spikes4 =c_spikes4';                        
                    end
                    id_cv_reps=0;
                    for idx=1:cv_reps_b
                            id_cv_reps=id_cv_reps+1;
                            % test and training set indices
                            c_spikes1_s=c_spikes1(randperm(nTrials),:);
                            c_spikes2_s=c_spikes2(randperm(nTrials),:); 

                            x_train  = [c_spikes1_s(3:end,:);c_spikes2_s(3:end,:)]; 
                            y_train = [ones(nTrials-2,1);ones(nTrials-2,1)*2];
                            x_test  = [c_spikes1_s(1:2,:);c_spikes2_s(1:2,:)];
                            y_test =[ones(2,1);ones(2,1)*2];
                            x_test_gen  = [c_spikes3;c_spikes4]; 
                            y_test_gen = [ones(size(c_spikes3,1),1);ones(size(c_spikes4,1),1)*2];
                            if flag_trial_shuffle==1
                                y_train=y_train(randperm(length(y_train)));
                            end
                            x_test(:,var(x_train)<0.000001)=[];  
                            x_test_gen(:,var(x_train)<0.000001)=[];  
                            x_train(:,var(x_train)<0.000001)=[];
               
                            %%% maltab classify
                            if ~isempty(x_train)
                                if training_method==2 %svm
                                    SVMModel=fitcsvm(x_train, y_train,'Standardize',true,'CacheSize','maximal','KernelScale','auto');
                                    [yHat,score] = predict(SVMModel,x_test);
                                    errs_generalize_real(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test(:) ~= yHat(:))/length(y_test);
                                    [yHat,score] = predict(SVMModel,x_test_gen);
                                    errs_generalize(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test_gen(:) ~= yHat(:))/length(y_test_gen);
                                elseif training_method==3 % ntL toolbox
                                    classifier=max_correlation_coefficient_CL;
                                    classifier = classifier.train(x_train', y_train);   
                                    [yHat, ~] = classifier.test(x_test');
                                    errs_generalize_real(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test(:) ~= yHat(:))/length(y_test);                                                                                                                    
                                    [yHat, decision_values] = classifier.test(x_test_gen');
                                    errs_generalize(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test_gen(:) ~= yHat(:))/length(y_test_gen);                                                                                
                                elseif training_method==1 %classify
                                    LDAModel=fitcdiscr(x_train, y_train,'SaveMemory','on');
                                    [yHat,~] = predict(LDAModel,x_test);                                
                                    errs_generalize_real(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test(:) ~= yHat(:))/length(y_test); 
                                    [yHat,~] = predict(LDAModel,x_test_gen);                                                                    
                                    errs_generalize(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test_gen(:) ~= yHat(:))/length(y_test_gen);                         
                                end   
                            end
                    end
                end
            end %%% end of generalize decoding pair
            %%  for non-generalize
            for c_decode_pair=1:N_decode_pair_b %%% for generalize
                for id_timebine=1:N_timebin_decoding   
                    c_spikes1=squeeze(c_spikes_side(:,decode_pair_b(c_decode_pair,1),id_timebine,:));
                    c_spikes2=squeeze(c_spikes_side(:,decode_pair_b(c_decode_pair,2),id_timebine,:));
                    c_spikes3=squeeze(c_spikes_side(:,decode_pair_b(c_decode_pair,3),id_timebine,:));
                    c_spikes4=squeeze(c_spikes_side(:,decode_pair_b(c_decode_pair,4),id_timebine,:));                   
                    if ~isvector(c_spikes1)
                        c_spikes1 =c_spikes1';                        
                        c_spikes2 =c_spikes2';
                        c_spikes3 =c_spikes3';
                        c_spikes4 =c_spikes4';                        
                    end
                    id_cv_reps=0;
                    for idx=1:cv_reps_b  
                        c_spikes1_s=c_spikes1(randperm(nTrials),:);
                        c_spikes2_s=c_spikes2(randperm(nTrials),:);
                        c_spikes3_s=c_spikes3(randperm(nTrials),:);
                        c_spikes4_s=c_spikes4(randperm(nTrials),:);
                        id_cv_reps=id_cv_reps+1;
                        x_train  = [c_spikes1_s(2:end,:);c_spikes2_s(2:end,:);c_spikes3_s(2:end,:);c_spikes4_s(2:end,:)]; 
                        y_train=[ones(nTrials-1,1);ones(nTrials-1,1)*2;ones(nTrials-1,1);ones(nTrials-1,1)*2];
                        x_test  = [c_spikes1_s(1,:);c_spikes2_s(1,:);c_spikes3_s(1,:);c_spikes4_s(1,:)];
                        y_test =[ones(nTest,1);ones(nTest,1)*2;ones(nTest,1);ones(nTest,1)*2];
                        if flag_trial_shuffle==1
                            y_train=y_train(randperm(length(y_train)));
                        end
                        x_test(:,var(x_train)<0.000001)=[];  
                        x_train(:,var(x_train)<0.000001)=[];
           
                        %%% maltab classify
                        if ~isempty(x_train)
                            if training_method==2 %svm
                                SVMModel=fitcsvm(x_train, y_train,'Standardize',true,'CacheSize','maximal','KernelScale','auto');
                                [yHat,score] = predict(SVMModel,x_test);
                            elseif training_method==3 % ntL toolbox
                                classifier=max_correlation_coefficient_CL;
                                classifier = classifier.train(x_train', y_train);   
                                [yHat, ~] = classifier.test(x_test');
                            elseif training_method==1 %classify
                                LDAModel=fitcdiscr(x_train, y_train,'SaveMemory','on');
                                [yHat,~] = predict(LDAModel,x_test);                                
                            end  
                            errs(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test(:) ~= yHat(:))/length(y_test);
                        end
                    end
                end
            end %%% non-end of generalize decoding pair
        end % end of loop for orientation
        errs_all{id_neurongroup}.errs=errs;
        errs_all{id_neurongroup}.errs_generalize_real=errs_generalize_real;
        errs_all{id_neurongroup}.errs_generalize=errs_generalize;
end
fprintf('Done.\n')

%%
save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_decording_error_',...
    nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},nameext_shuffle{flag_trial_shuffle},num2str(decoding_timebin),'.mat']),'errs_all','-v7.3');   

%%
% figure;
% % aa=squeeze(mean(errs_generalize(4,:, :,1,:),3));
% % aa=squeeze(mean(errs(1,:, :,1,:),3));
% image(bincenter_psth_downsample(1:40),1:20,1-aa','CDataMapping','scaled')
% hold on
% colorbar
%         ax=gca;
%         ax.CLim=([0.4,0.9]);
%                 ax.YDir='Normal';
%             xLimits = get(gca,'XLim');
% %             if flag_exis_layer==1
% %                 for temp_i=1:length(cell_layer_idx_border_included)
% %                     plot(xLimits,[cell_layer_idx_border_included(temp_i)+0.5,cell_layer_idx_border_included(temp_i)+0.5],'k--','LineWidth',layerline_width)
% %                     hold on
% %                 end 
% %             end
% % plot(bincenter_psth_downsample,1-aa(:,4),'r')
% hold on

% % 
% SVMModel=fitcsvm(x_train, y_train,'Standardize',true,'CacheSize','maximal','Solver','L1QP','OptimizeHyperparameters','auto',...
%    'HyperparameterOptimizationOptions',struct('AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0));
