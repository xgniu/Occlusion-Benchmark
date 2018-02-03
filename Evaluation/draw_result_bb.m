close all;clear;clc;warning off all;

addpath('.\util');
pathRes = '.\results\';% The folder containing the tracking results
pathDraw = '.\tmp\';% The folder that will stores the images with overlaid bounding box

seqs = config_occ_seqs();
trks = config_trackers();

LineWidth = 2;
plotSetting();

for index_seq=1:length(seqs)
    seq = seqs{index_seq};
    seq_name = seq.name;
    seq_length = seq.endFrame-seq.startFrame+1; %size(rect_anno,1);
    
    resultsAll=[];
    trackerNames=[];
    for index_algrm = 1:length(trks)
        algrm = trks{index_algrm};
        name = algrm.name;
        trackerNames{index_algrm}=name;
               
        fileName = [pathRes seq_name '_' name '.mat'];    
        results = load(fileName);        
        res = results.bboxes;
        resultsAll{index_algrm} = res;
    end
    
    d = dir(['D:\Visual Tracking\Benchmark\OCC\',seq.name,'\*.jpg']);
    for i = 1 : size(d, 1)
        img{i} = [d(i).folder, '\', d(i).name];
    end
        
    pathSave = [pathDraw seq_name '\'];
    if ~exist(pathSave,'dir')
        mkdir(pathSave);
    end
    
    for i=1:seq_length
        im = imread(img{i}); 
        imshow(im);
        text(10, 15, ['#' num2str(i)], 'Color','y', 'FontWeight','bold', 'FontSize',24);        
        for j=1:length(trks)
            LineStyle = plotDrawStyle{j}.lineStyle;
            rectangle('Position', resultsAll{j}(i,:), 'EdgeColor', plotDrawStyle{j}.color, 'LineWidth', LineWidth,'LineStyle',LineStyle); 
        end        
        imwrite(frame2im(getframe(gcf)), [pathSave num2str(i) '.png']);
    end
    clf
end
