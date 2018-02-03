clear;close all;addpath('.\util');

fig_path = '.\figs\';
perf_mat_path = '.\perfMat\';

plot_draw_style = {
    struct('color',[1,0,0],'lineStyle','-'), struct('color',[0,1,0],'lineStyle','--'), struct('color',[0,0,1],'lineStyle',':'),...
    struct('color',[0,0,0],'lineStyle','-'), struct('color',[1,0,1],'lineStyle','--'), struct('color',[0,1,1],'lineStyle',':'),...
    struct('color',[0.5,0.5,0.5],'lineStyle','-'), struct('color',[136,0,21]/255,'lineStyle','--'),...%dark red
    struct('color',[255,127,39]/255,'lineStyle',':'), struct('color',[0,162,232]/255,'lineStyle','-'),...%Turquoise
    };

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

gt_all = [];
for idxSeq = 1:num_seqs
    seq = seqs{idxSeq};
    gt = dlmread(['D:\Visual Tracking\Benchmark\OCC\' seq.name '\groundtruth_rect.txt']);
    gt_ = zeros(size(gt,1),4);
    for i=1:size(gt,1)
        region = gt(i,:); %initial rectangle [x,y,width, height]
        if numel(region) > 4
            x1 = round(min(region(1:2:end)));
            x2 = round(max(region(1:2:end)));
            y1 = round(min(region(2:2:end)));
            y2 = round(max(region(2:2:end)));
            region = round([x1, y1, x2 - x1, y2 - y1]);
        else
            region = round(region);
        end        
        gt_(i,:) = region;
    end
    gt_all = [ gt_all; gt_];
end

%% Success Plot:
figure;
legends = [];
for idxTrk = 1:num_trackers
    tracker = trackers{idxTrk};
    res_all = [];
    for idxSeq = 1:num_seqs
        seq = seqs{idxSeq};
        results = load(['.\results\' seq.name '_' tracker.name '.mat']);% load seq_tracker_results
        res_all = [ res_all; results.bboxes];
    end
    
    err_overlap = calc_overlap_err(res_all, gt_all);
    
    overlap_thresholds = 0:0.05:1;
    success_scores = zeros(1,length(overlap_thresholds));   
    for tIdx=1:length(overlap_thresholds)
        success_scores(1,tIdx) = sum( err_overlap > overlap_thresholds(tIdx) ) / size(res_all,1);
    end
    
    hold on;
    plot(overlap_thresholds,success_scores,'color',plot_draw_style{idxTrk}.color, 'lineStyle', ...
            plot_draw_style{idxTrk}.lineStyle,'lineWidth', 2);
    legends{idxTrk} = [tracker.name ' [' sprintf('%.3f', success_scores(11)) ']'];
end
legend(legends,'fontsize',10);


%% Precision Plot:
figure;
legends = [];
for idxTrk = 1:num_trackers
    tracker = trackers{idxTrk};
    res_all = [];
    for idxSeq = 1:num_seqs
        seq = seqs{idxSeq};
        results = load(['.\results\' seq.name '_' tracker.name '.mat']);% load seq_tracker_results
        res_all = [ res_all; results.bboxes];
    end
    
    ncle = calc_center_err(res_all, gt_all);
    
    ncle_thresholds = 0:0.05:1;
    precision_scores = zeros(1,length(ncle_thresholds));   
    for tIdx=1:length(ncle_thresholds)
        precision_scores(1,tIdx) = sum( ncle < ncle_thresholds(tIdx) ) / size(res_all,1);
    end
    
    hold on;
    plot(ncle_thresholds,precision_scores,'color',plot_draw_style{idxTrk}.color, 'lineStyle', ...
            plot_draw_style{idxTrk}.lineStyle,'lineWidth', 2);
    legends{idxTrk} = [tracker.name ' [' sprintf('%.3f', precision_scores(11)) ']'];
end
legend(legends,'fontsize',10);