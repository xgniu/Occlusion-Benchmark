function err_center = calc_center_err(bboxes, gt)
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
    
    center = [ bboxes(:,1)+(bboxes(:,3)-1)/2 bboxes(:,2)+(bboxes(:,4)-1)/2 ];
    center_gt = [ gt_(:,1)+(gt_(:,3)-1)/2 gt_(:,2)+(gt_(:,4)-1)/2];
    d1 = abs(center(1:seq_length,1) - center_gt(1:seq_length,1)) ./ gt(1:seq_length,3);
    d2 = abs(center(1:seq_length,2) - center_gt(1:seq_length,2)) ./ gt(1:seq_length,4);
    err_center = min( max(d1,d2), 1 );
    err_center(~idx) = -1;
    
end
