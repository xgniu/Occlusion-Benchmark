clear
close all;

addpath('.\util');

fig_path = '.\figs\';
perf_mat_path = '.\perfMat\';

plot_draw_style = {
    struct('color',[1,0,0],'lineStyle','-'), struct('color',[0,1,0],'lineStyle','--'), struct('color',[0,0,1],'lineStyle',':'),...
    struct('color',[0,0,0],'lineStyle','-'), struct('color',[1,0,1],'lineStyle','--'), struct('color',[0,1,1],'lineStyle',':'),...
    struct('color',[0.5,0.5,0.5],'lineStyle','-'), struct('color',[136,0,21]/255,'lineStyle','--'),...%dark red
    struct('color',[255,127,39]/255,'lineStyle',':'), struct('color',[0,162,232]/255,'lineStyle','-'),...%Turquoise
    }; % 1x10 cell array
rankNum = 10; % show top-10 rank trackers

%% Load infos of seqs and trackers
seqs = config_occ_seqs();
num_seqs = length(seqs);
name_seqs = cell(num_seqs,1); % name of all seqs
for idxSeq=1:num_seqs
    t = seqs{idxSeq};
    name_seqs{idxSeq}=t.name;
end

trackers = config_trackers();
num_trackers = length(trackers);
name_trackers = cell(num_trackers,1); % name of all trackers
for idxTrk=1:num_trackers
    t = trackers{idxTrk};
    name_trackers{idxTrk}=t.namePaper;
end

%% 
eval_types = {'OPE'};

%% Success Plot:
metric_type='Success';
thresholds = 0:0.05:1;
rankIdx = 11; % threshold(rankIdx)=0.500, 50% overlap ratio, for ranking
xLabelName = 'Overlap threshold';
yLabelName = 'Success rate';
rankingType = 'AUC';

for j = 1:length(eval_types) % for OPE,SRE,TRE:
    eval_type = eval_types{j}; %SRE, TRE, OPE
    plot_type = [metric_type '_' eval_type];
    
    gen_success_mat(seqs, trackers, eval_type, name_trackers, perf_mat_path);   
    
    dataName = [perf_mat_path 'SuccessMat_' num2str(num_trackers) 'algs_' eval_type '.mat'];
    success_mat = load(dataName);
    num_trackers = size(success_mat.avg_success,1);    
    if rankNum > num_trackers || rankNum <0 % if numTrk is less than rankNum(10),
        rankNum = num_trackers;
    end
    
    titleName = ['Success Plots of ' eval_type];
    figName= [fig_path plot_type '_' rankingType];
    % draw and save the overall performance plot
    idxSeqSet = 1:length(seqs);
    plot_draw_save(num_trackers,plot_draw_style,success_mat.avg_success,idxSeqSet,rankNum,rankingType,rankIdx,...
        name_trackers,thresholds,titleName, xLabelName,yLabelName,figName);

end


%% Precision Plot
metric_type = 'Precision';

thresholds = 0:0.05:1;
rankIdx = 11; % thresholdSet(rankIdx)=20, 20pixel location error, for ranking
xLabelName = 'Location error threshold';
yLabelName = 'Precision';
rankingType = 'threshold';

for j=1:length(eval_types) % for OPE,SRE,TRE:
    eval_type = eval_types{j};%SRE, TRE, OPE
    plot_type = [metric_type '_' eval_type];
    
    gen_precision_mat(seqs, trackers, eval_type, name_trackers, perf_mat_path);
      
    dataName=[perf_mat_path 'PrecisionMat_' num2str(num_trackers) 'algs_' eval_type '.mat'];
    precision_mat = load(dataName);
    num_trackers = size(precision_mat.avg_precision,1);
    
    if rankNum > num_trackers || rankNum <0 % if numTrk is less than rankNum(10),
        rankNum = num_trackers;
    end
    
    titleName = ['Precision plots of ' eval_type];
    figName= [fig_path plot_type '_' rankingType];
    % draw and save the overall performance plot
    idxSeqSet = 1:length(seqs);    
    plot_draw_save(num_trackers,plot_draw_style,precision_mat.avg_precision,idxSeqSet,rankNum,rankingType,rankIdx,...
        name_trackers,thresholds,titleName, xLabelName,yLabelName,figName);

end
