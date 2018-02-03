function err_overlap = calc_overlap_err(bboxes, gt)
    % bboxes: A run of the tracker on the (sub)seq.
    % gt: corresponding groundtruth
    % errCoverage: overlap ratio error, for Success Plot, (seq_length,1)
    % errCenter: center location error, for Precision Plot, (seq_length,1)
    
    % only work for 'rect' types of results
    
    %% load results and validate
    seq_length = min(size(bboxes,1),size(gt,1));
    gt_ = zeros(seq_length,4);
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
    
    bboxes(1,:) = gt_(1,:);
   
    index = gt_>0;
    idx=(sum(index,2)==4); % idx of valid rows of gt
    idx = idx(1:seq_length);

    tmp = calc_rect_int(bboxes(idx,:),gt_(idx,:));    
    err_overlap = -ones(length(idx),1);
    err_overlap(idx) = tmp;
    
end
