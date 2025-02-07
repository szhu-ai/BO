% This is the first step to analyze neuropixel recording data
% (1) load event/trigger channels, get event time based on raw voltage,
% mannully modify based on the number of events(stim) you already have
% (2) calculate spike counts/ FR, save it into matrix for later.

% Base on example script for some of the functions in the spikes repository. 
%
% See https://github.com/cortex-lab/neuropixels/wiki/Other_analysis_methods
% for more explanation. 
clear all
clc
%% add the repositories to your path
srcPath ='G:\npix\code\xm_kilosort';
addpath(genpath(fullfile(srcPath,'KiloSort'))); % path to kilosort folder
addpath(genpath(fullfile(srcPath,'npy-matlab'))); % path to npy-matlab scripts
addpath(genpath(fullfile(srcPath,'spikes')));
% addpath(genpath(fullfile(srcPath,'tantalus')));
addpath(genpath(fullfile(srcPath,'export_fig')));
% addpath(genpath(fullfile(srcPath,'workups/runKilosort')));
codePath ='G:\npix\code';
addpath(genpath(fullfile(codePath,'decoding_glmnet')));
%% set paths for where to find your data, please edit this part accordingly
rootpath='G:\npix\';
testtype='bo';
recordingDate='20190305'; %%% for example, 20181206
recordingSession='R02';

bankid=1;
event_on_time_exist=1;  %%% set to zero if you haven't got the event/trigger time

myKsDir = 'G:\Quicky_kilo2_bo\';
myresultPath=[rootpath,'spike\BO\'];
myresultPath_psth=[rootpath,'result\BO\'];
myEventDir=[rootpath,'event_time\bo\']; %%% location for event on time 
myStimPath=[rootpath,'event\'];
%% load condition id for all the trials
if strcmp(testtype,'bo') && bankid==2
    load([myStimPath,recordingDate,'\',recordingSession,'\','bo_all2.mat']);   %% for data 181206, we have two orientation test for two banks
    load([myStimPath,recordingDate,'\',recordingSession,'\','bo_trialmat2.mat']);
elseif strcmp(testtype,'bo') && bankid==1
    load([myStimPath,recordingDate,'\',recordingSession,'\','bo_all.mat']);   %% for data 181206, we have two orientation test for two banks
    load([myStimPath,recordingDate,'\',recordingSession,'\','bo_trialmat.mat']);
%     load([myStimPath,recordingDate,'\',recordingSession2,'\','bo_trialmat.mat']);
end
c_eye=5;
c_ori=6;
c_side=7;
c_lc=8;
c_sz=9;
trialmat(trialmat(:,c_side)==2,c_lc)=3-trialmat(trialmat(:,c_side)==2,c_lc);

%% set parameters for the timing counting window for different analysis
spikecount_window=[-100,500]; %% in milliseconds, relative to stim onset(0)
time_bin_sz=50;
time_bin=spikecount_window(1):time_bin_sz:spikecount_window(2);
time_bin_center=(time_bin(1)+time_bin_sz/2):time_bin_sz:(time_bin(end)-time_bin_sz/2);
N_timebin=length(time_bin)-1;

load([myresultPath,recordingDate,'_',recordingSession,'\',testtype,num2str(bankid),'_sp.mat']);
load([myresultPath,recordingDate,'_',recordingSession,'\',testtype,num2str(bankid),'_Cluster.mat']);
load([myresultPath,recordingDate,'_',recordingSession,'\',testtype,num2str(bankid),'_spikes.mat']);

%% compute the firing rate for each condition
% S_event= load([myEventDir,recordingDate,'\',recordingSession,'\','EventOn_',testtype,num2str(bankid),'.mat']); % a vector of times in seconds of some event to align to 
% event_on_time=S_event.event_on_time;
% N_cluster=length(cluster.depthsorted_id);
% FR_all=zeros(N_stim,N_timebin,N_cluster);
% 
% for i=1:length(cluster.depthsorted_id)
%     if cluster.depthsorted_label(i)~=3
%         sptime=sp.st(sp.clu==cluster.depthsorted_id(i));
%         S1=arrayfun(@(x) histcounts(sptime*1000-x,time_bin), event_on_time*1000,'UniformOutput',false);
%         FR_all(:,:,i)=1000*cell2mat(S1)./time_bin_sz;
%     end
% end 
% 
% 
% idx=0;
% spikes=zeros(N_cluster,N_condition,N_timebin,stim_repetitions);
% for ori_id=1:N_ori
%     for side_id=1:N_side
%         for lc_id=1:N_lc
%             for sz_id=1:N_sz
%                 idx=idx+1;
%                 FR_currentcond=FR_all(trialmat(:,c_ori)==ori_id&trialmat(:,c_side)==side_id&trialmat(:,c_lc)==lc_id&trialmat(:,c_sz)==sz_id,:,:);
%                 spikes(:,idx,:,:)=permute(FR_currentcond,[3,2,1]);
%             end
%         end
%     end    
% end
% save([myresultPath,recordingDate,'_',recordingSession,'\',testtype,num2str(bankid),'_spikes.mat'],'spikes');
%% 

current_ori=4; %%% 1,2,3,4
%% Select a dataset on which to perform the analysis 
times = time_bin_center;
spk = spikes;
a=[1 1 1 1 2 2 2 2];
conditions.side=repmat(a',[4,1]);
b=[1 2 3 4 1 2 3 4];
conditions.contrast=repmat(b',[4,1]);
c=[1 2 3 4];
conditions.ori=sort(repmat(c',[8,1]));
% use only active cells, threshold .1 Hz across all stimuli
active = mean(sum(mean(spk(:,:,times>=0 & times<=500,:),4),3)/0.5,2) > 1;
spk = spk(active,:,:,:);

spk=spk(:,:,:,:);
layer='all';



side =[1 2];
contrast = [2 4]; nCon = length(contrast);

%% Compute features: Filtering with a 50 ms box function
winSize = 100;
binWidth = 20;
n = fix(winSize/binWidth/2);
win = ones(2*n+1,1);                    % construct window
win = reshape(win,1,[])/length(win);        

fprintf('Feature extration.\n')
for n = 1:length(contrast)
  for k = 1:length(side)
    % find the condition
    idx = find(cat(1,conditions.side)==side(k) & ...
                cat(1,conditions.contrast)==contrast(n)&conditions.ori==current_ori);
    
    feat{n,k} = filter(win,1,permute(spk(:,idx,:,:),[1 3 4 2]),[],2); %#ok<SAGROW>

  end
end
[nNeurons nBins nTrials] = size(feat{1,1});

%% Run cross-validated decoding on two orientations for both contrasts
 
numClasses = length(side);
yTrain = kron([1:numClasses]', ones(nTrials,1)); % assume each class has the same number of trials, # x 1
  
nTrials = length(yTrain); % overwrite previous one with total trial numbers across all the classes

cv_reps = 20;                  % number of repetitions for each cv
                                % decrease for better performance
nTest = ceil(nTrials * .2);        % 20% of the trials for testing, 
                                % 80% for training during cross-validation
% ori1 = 1; ori2 = 3;             % just choosing two orientations for speed

% nTrain = nTrials - nTest; % cv training trials, used for performance measure

% Elastic net mixing parameter
alpha = 0.8; % 0.8, adjustable [0,1], alpha=0, L2; alpha=1, lasso

% Number of CV folds in validating regularization parameter lambda
nfolds = 5;

% Calculate regularization path
% verbose = false;
options = glmnetSet();
options.alpha = alpha;
        
% err = NaN(nCon,nBins);
fprintf('Classification.\n')
% warning off         %#ok<WNOFF>     % glmnet toolbox throws lots of warnings 
for n = 1:nCon
    fprintf('Contrast %d\n',contrast(n))
    for b=1:nBins
        if mod(b,5)==0,fprintf('.'),end        
 
        % Data preparation
        ddd=[];
        for k = 1:length(side),
            dtmp = feat{n,k}(:,b,:);
            ddd = cat(3,ddd, dtmp);
        end
        DAT = permute(ddd, [1 3 2])'; % # x nNeurons     
        
        nTrials = size(DAT, 1);
        
        % generate partitions for cross-validation
        warning off
        cvt1 = getCvTrials(nTrials,nTest,cv_reps);
        warning on
        
        for c=1:cv_reps
            % test and training set indices
            idx_test1 = cvt1(c,:);     
            idx_train1 = setdiff(1:nTrials,idx_test1);
            
            % prepare training and test set        
            x_test  = DAT(idx_test1,:);
            y_test = yTrain(idx_test1);
            x_train = DAT(idx_train1,:);
            y_train = yTrain(idx_train1);
            
            warning off %#ok
            % cvErr = cvglmnetMulticlass2(x_train, y_train, nfolds, [], 'class','multinomial', options, verbose);
            cvErr = cvglmnet(x_train, y_train,'multinomial',options,'class',nfolds);
            % cvglmnetPlot(cvErr)
            warning on %#ok

            % Choose model with lowest CV error
            % lambda = cvErr.lambda_min;
            lambda = cvErr.lambda_1se;
            best = find(cvErr.lambda == lambda);
            beta = [cvErr.glmnet_fit.a0(:,best), zeros(numClasses, nNeurons)];
            for k = 1 : numClasses
                beta(k,2:end) = cvErr.glmnet_fit.beta{k}(:,best)';
            end
            
            % testing decoder
             yHat = repmat(beta(:,1), 1, nTest) + beta(:,2:end)*x_test';
            [~,yHat] = max(yHat);

            errs(n,b, c) = sum(y_test(:) ~= yHat(:))/nTest;           
        end % c        
        
    end
    fprintf('\n')
end % n
% warning on %#ok<WNON>

fprintf('Done.\n')

err= mean(errs, 3);

%% Plot fraction correct classification performance as function of time 

figure
% set(gcf,'position',[100 100 400 400])
plot(times,1-err(1,:),'color',[1 0 0])
hold on
plot(times,1-err(2,:),'color',[0 0 1])
hold on
axis square

set(gca,'box','off')
xlabel('Time (ms)')
ylabel('Percent correct')
ca = axis;
% axis([-200 700 0. ca(4)])
axis([-200 500 0 1])

line([-200 500],[1/numClasses 1/numClasses],'color','k','linestyle',':')
% line([0 0],[0. ca(4)],'color','k','linestyle',':')
line([0 0],[0 1],'color','k','linestyle',':')


l = legend({'LC0,SZ8','LC1,SZ8'});
set(l,'box','off')
ax=gca;
ax.XLabel.FontWeight='Bold';
ax.YLabel.FontWeight='Bold';
filename=[recordingDate,'_',recordingSession,'_decoding_',num2str(current_ori),'_',layer];
print('-dpng', [myresultPath_psth,filename,'.png'], '-r300');

%%
