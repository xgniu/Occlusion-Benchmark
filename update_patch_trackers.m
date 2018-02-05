function occ = update_patch_trackers(track_output,im)
    % track_output: (x,y,width,height)    
    global patches;
    num_patches = numel(patches);
    global occ_bboxes;
    
    track_output_center = round( track_output([1 2])+track_output([3 4])/2 );%(x,y)
    
    % im: the current frame
    if size(im,3)~=1
        imgray=single(rgb2gray(im));
    else
        imgray=single(im);
    end
    
    vs = track_output(2):track_output(2)+track_output(4); %vertical scale
    hs = track_output(1):track_output(1)+track_output(3); %horizontal scale
    vs(vs < 1) = 1;
    hs(hs < 1) = 1;
    vs(vs > size(im,1)) = size(im,1); % height
    hs(hs > size(im,2)) = size(im,2); % width
    vs = round(vs);
    hs = round(hs);
    mask = zeros(size(im,1),size(im,2));
    mask(vs,hs) = 1;
    total_area = sum(mask(:));
    
    for i = 1:num_patches
        patch = patches{i};
        [patch.box,response] = patch.tracker.track(imgray);
        psr = PSR(response);
        overlap_ratio = rectint(patch.box,track_output) / (patch.box(3)*patch.box(4));
        
        if psr>50 && overlap_ratio>0
            patch_box = patch.box;
            vs = patch_box(2):patch_box(2)+patch_box(4); %vertical scale
            hs = patch_box(1):patch_box(1)+patch_box(3); %horizontal scale
            vs(vs < 1) = 1;
            hs(hs < 1) = 1;
            vs(vs > size(im,1)) = size(im,1);
            hs(hs > size(im,2)) = size(im,2);
            vs = round(vs);
            hs = round(hs);
            mask(vs,hs) = 0;            
        end
        
        if psr > 50 && overlap_ratio > 0.3
            occ_bboxes(i) = 1;
            patch.tracker.update();
        else
            patch = reset_patch(imgray,track_output_center([2 1]),track_output([4 3]),patch.tag);
        end
        patches{i} = patch;
        
    end
    
    remain_area = sum(mask(:));
    occ = (remain_area/total_area) < 0.85;

end