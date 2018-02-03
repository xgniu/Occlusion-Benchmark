function [img, init_center, target_sz, gt] = load_seq_info(base_path,seqName)
    % Load video info.
    
    data_path=[base_path seqName];
    
    groundtruth_rect = importdata([data_path '/groundtruth_rect.txt']);
    gt = zeros(size(groundtruth_rect,1),4);  %[x,y,width, height]
    
    for i=1:size(groundtruth_rect,1)
        region = groundtruth_rect(i,:); %initial rectangle [x,y,width, height]

        if numel(region) > 4
            x1 = round(min(region(1:2:end)));
            x2 = round(max(region(1:2:end)));
            y1 = round(min(region(2:2:end)));
            y2 = round(max(region(2:2:end)));
            region = round([x1, y1, x2 - x1, y2 - y1]);
        else
            region = round(region);
        end
        
        gt(i,:) = region;
    end
    
    initstate = gt(1,:); %initial rectangle [x,y,width, height]
    %center of the target, (row, col)
    init_center = [initstate(2)+initstate(4)/2 initstate(1)+initstate(3)/2];
    %initial size of the target, (rows, cols)
    target_sz = [initstate(4) initstate(3)];
    
    d = dir([data_path, '/*.jpg']);
    for i = 1 : size(d, 1)
        img{i} = [d(i).folder, '\', d(i).name];
    end
    
end