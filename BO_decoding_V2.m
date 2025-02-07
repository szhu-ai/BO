
%% set parameters for the timing counting window for different analysis
flag_trial_shuffle=2;
nameext_shuffle={'shuffle_','real_'};

flag_singleneuron=2;
nameext_neuron={'single_','population_'};

flag_singletime=2;
nameext_time={'onetime','alltime'};  

training_method=1; 
nameext_method={'classify_','svm_'};  

%% compute or load the firing rate for each condition
[~,sortidx]=sort(id_condition);
aa=reshape(data(sortidx,:,:),N_repetition,N_condition,N_Bins,N_cluster);
spikes_all_temp=permute(aa,[4,2,3,1]); %%(N_cluster,N_condition,N_timebin,N_repetition);
clear aa
%%    
window_decoding=[0,0.35];
window_decoding_idx=find(binborder_psth>=window_decoding(1),1,'first') :1: find(binborder_psth<=window_decoding(2),1,'last')-1;
bincenter_decoding=bincenter_psth(window_decoding_idx);
x_downsample=1:10:length(window_decoding_idx);

spikes_decoding=spikes_all_temp(Includedidx_all,:,window_decoding_idx,:)./psth_binSize;
if flag_singletime==1
    spikes_downsample=mean(spikes_decoding,3);
    N_timebin_decoding=1;
else  % smooth data within 100 ms and then downsample by 20 times
    spikes_sm=smoothdata(spikes_decoding,3,'gaussian',100);
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
decode_pair=[1 2; 3 4; 1 3; 2 4];
decode_pair2=[1,2,3,4; 3,4,1,2;1,3,2 4; 2,4,1,3];
[N_decode_pair,numClasses]=size(decode_pair);
N_decode_pair2=size(decode_pair2,1);
nTest=1;                             % leave one out %   nTest = ceil(nTrials * .2);  
nTrials=size(spikes_downsample,4);  %% 2 binary class compariasion
cvt1 = nchoosek(1:nTrials,nTest);  
ncvt1=size(cvt1,1);
cv_reps=min(ncvt1,1000); 
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

        errs=zeros(N_decode_pair,N_timebin_decoding,cv_reps,n_ori,'single')+0.5;
        errs_generalize=zeros(N_decode_pair2,N_timebin_decoding,nTrials*2,n_ori,'single')+0.5;

        fprintf('Group %d \n',id_neurongroup)
        for id_ori = 1:4
            if flag_singleneuron==1
                c_spikes_side=spikes_downsample(id_neurongroup,(id_ori-1)*8+[2,4,6,8],:,:);
            elseif flag_singleneuron==2
                c0=id_neuron_downsample(id_neurongroup)-4;
                c1=id_neuron_downsample(id_neurongroup)+4;
                c0=max(c0,1);
                c1=min(c1,N_neuron_group);
                c_spikes_side=spikes_downsample(c0:c1,(id_ori-1)*8+[2,4,6,8],:,:);
            end
            %%  for generalize
            for c_decode_pair=1:4 %%% for generalize
                c_spikes_train=c_spikes_side(:,decode_pair2(c_decode_pair,1:2),:,:);
                c_spikes_test=c_spikes_side(:,decode_pair2(c_decode_pair,3:4),:,:);
                x_train_mean=mean(c_spikes_train,[2,3,4]);
                x_train_var= var(c_spikes_train,0,[2,3,4]);                
                x_train_mean=x_train_mean';
                x_train_var=x_train_var';
                x_train_var(x_train_var==0)=1;

                for id_timebine=1:N_timebin_decoding   
                    c_spikes1=squeeze(c_spikes_side(:,decode_pair2(c_decode_pair,1),id_timebine,:));
                    c_spikes2=squeeze(c_spikes_side(:,decode_pair2(c_decode_pair,2),id_timebine,:));
                    c_spikes3=squeeze(c_spikes_side(:,decode_pair2(c_decode_pair,3),id_timebine,:));
                    c_spikes4=squeeze(c_spikes_side(:,decode_pair2(c_decode_pair,4),id_timebine,:));                   
                    if ~isvector(c_spikes1)
                        c_spikes1 =c_spikes1';                        
                        c_spikes2 =c_spikes2';
                        c_spikes3 =c_spikes3';
                        c_spikes4 =c_spikes4';                        
                    end
                    x_train  = [c_spikes1;c_spikes2];
                    x_train_zscore=(x_train-repmat(x_train_mean,size(x_train,1),1))./repmat(x_train_var,size(x_train,1),1);
                    for j=1:size(x_train_zscore,2)
                        x_train_zscore(abs(x_train_zscore(:,j))>=4,j)=mean(x_train_zscore(abs(x_train_zscore(:,j))<4,j));
                    end
                    y_train = [ones(nTrials,1);ones(nTrials,1)*2];
                    if flag_trial_shuffle==1
                        y_train1 = [ones(nTrials/2,1);ones(nTrials/2,1)*2];
                        y_train=[y_train1(randperm(nTrials));y_train1(randperm(nTrials))];
                    end

                    x_test = [c_spikes3;c_spikes4];
                    x_test_zscore=(x_test-repmat(x_train_mean,size(x_test,1),1))./repmat(x_train_var,size(x_test,1),1);
                    y_test = [ones(nTrials,1);ones(nTrials,1)*2];

                    x_test_zscore(:,var(x_train_zscore)<0.001)=[];  
                    x_train_zscore(:,var(x_train_zscore)<0.001)=[];

                    %%% maltab classify
                    if ~isempty(x_train_zscore)
%                         if training_method==2 %svm
%                             SVMModel=fitcsvm(x_train_zscore, y_train,'Standardize',true,'CacheSize','maximal','KernelScale','auto');
%                             [yHat,score] = predict(SVMModel,x_test_zscore);
%                             errs_generalize(c_decode_pair,id_timebine, :,id_ori)=y_test(:) ~= yHat(:);
%                         else
                            LDAModel=fitcdiscr(x_train_zscore, y_train,'SaveMemory','on');
                            [yHat,~] = predict(LDAModel,x_test_zscore);                                
                            errs_generalize(c_decode_pair,id_timebine, :,id_ori)=y_test(:) ~= yHat(:);                         
%                         end   
                    end
                end
            end %%% end of generalize decoding pair



            % for non-generalize
            for c_decode_pair=1:N_decode_pair
                for id_timebine=1:N_timebin_decoding                   
                    c_spikes1=squeeze(c_spikes_side(:,decode_pair(c_decode_pair,1),id_timebine,:));
                    if ~isvector(c_spikes1)
                        c_spikes1 =c_spikes1';
                    end
                    c_spikes2=squeeze(c_spikes_side(:,decode_pair(c_decode_pair,2),id_timebine,:));
                    if ~isvector(c_spikes2)
                        c_spikes2 =c_spikes2';
                    end   

                    % generate partitions for cross-validation
                    id_cv_reps=0;
                    for idx_trial_class1=1:ncvt1
%                         for idx_trial_class2=1:ncvt1
                            id_cv_reps=id_cv_reps+1;
                            % test and training set indices
                            idx_test_class1 = cvt1(idx_trial_class1,:);     
                            idx_train_class1 = setdiff(1:nTrials,idx_test_class1);     
                            idx_test_class2 = cvt1(idx_trial_class1,:);     
                            idx_train_class2 = setdiff(1:nTrials,idx_test_class2);                               
                            % prepare training and test set        
                            x_test  = [c_spikes1(idx_test_class1,:);c_spikes2(idx_test_class2,:)];
                            x_train  = [c_spikes1(idx_train_class1,:);c_spikes2(idx_train_class2,:)];
                            x_train_zscore=zscore(x_train);
                            for j=1:size(x_train,2)
                                x_train(abs(x_train_zscore(:,j))>=4,j)=mean(x_train(abs(x_train_zscore(:,j))<4,j));
                            end
                            x_train_mean=mean(x_train);
                            x_train_var= var(x_train);
                            x_train_var(x_train_var==0)=1;
                            x_train_zscore=zscore(x_train);
                            x_test_zscore=(x_test-repmat(x_train_mean,size(x_test,1),1))./repmat(x_train_var,size(x_test,1),1);
        
                            x_train_var=var(x_train);
                            x_train_zscore(:,x_train_var<0.1)=[];
                            x_test_zscore(:,x_train_var<0.1)=[];              

                            y_test =[ones(nTest,1);ones(nTest,1)*2];
                            y_train = [ones(length(idx_train_class1),1);ones(length(idx_train_class2),1)*2];
                            if flag_trial_shuffle==1
                                y_train=y_train(randperm(size(y_train,1)),1);                   
                            end
                            %% maltab classify
                            if ~isempty(x_train_zscore)
%                                 if training_method==2 %svm
%                                     SVMModel=fitcsvm(x_train_zscore, y_train,'Standardize',true,'CacheSize','maximal','KernelScale','auto');
%                                     [yHat,score] = predict(SVMModel,x_test_zscore);
%                                     errs(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test(:) ~= yHat(:))/length(y_test);
%                                 else
                                    LDAModel=fitcdiscr(x_train_zscore, y_train,'SaveMemory','on');
                                    [yHat,~] = predict(LDAModel,x_test_zscore);                                
                                    errs(c_decode_pair,id_timebine, id_cv_reps,id_ori)=sum(y_test(:) ~= yHat(:))/length(y_test);                                                                                
%                                 end
                            end
%                         end
                    end       
                end %%% end of for loop timebin
            end %%% end of non-generalize decoding pair
        end % end of loop for orientation
%         errs_all{id_neurongroup}.errs=errs;
        errs_all{id_neurongroup}.errs_generalize=errs_generalize;
end
fprintf('Done.\n')

%%
save(fullfile(myresultPath,recordingDate,recordingSession,[testtype,num2str(bankid),'_decording_error_',...
    nameext_method{training_method},nameext_neuron{flag_singleneuron},nameext_time{flag_singletime},nameext_shuffle{flag_trial_shuffle},'generalize.mat']),'errs_all','-v7.3');   

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
