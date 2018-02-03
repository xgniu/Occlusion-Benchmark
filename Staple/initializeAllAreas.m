function seq = initializeAllAreas(im, params, seq)

	% we want a regular frame surrounding the object
	avg_dim = sum(seq.init_target_sz)/2;
	% size from which we extract features for correlation filter
	bg_area = round(seq.init_target_sz + avg_dim);
	% pick a "safe" region smaller than bbox to avoid mislabeling
	fg_area = round(seq.init_target_sz - avg_dim * params.inner_padding);
	% saturate to image size
	if(bg_area(2)>size(im,2)), bg_area(2)=size(im,2)-1; end
	if(bg_area(1)>size(im,1)), bg_area(1)=size(im,1)-1; end
	% make sure the differences are a multiple of 2 (makes things easier later in color histograms)
    % make the bg_area has the same parity as target_sz
	seq.bg_area = bg_area - mod(bg_area - seq.init_target_sz, 2);
    % make the fg_area has the same parity as bg_area
	seq.fg_area = fg_area + mod(seq.bg_area - fg_area, 2);

	% Compute the rectangle with (or close to) params.fixedArea and
	% same aspect ratio as the target bbox
	seq.area_resize_factor = sqrt(params.fixed_area/prod(seq.bg_area));
	seq.norm_bg_area = round(seq.bg_area * seq.area_resize_factor);
	
    % Correlation Filter (HOG) feature space
	% It smaller that the norm bg area if HOG cell size is > 1
	seq.cf_response_size = floor(seq.norm_bg_area / params.hog_cell_size);
	
    % given the norm BG area, which is the corresponding target w and h?
 	norm_target_sz_w = 0.75*seq.norm_bg_area(2) - 0.25*seq.norm_bg_area(1);
 	norm_target_sz_h = 0.75*seq.norm_bg_area(1) - 0.25*seq.norm_bg_area(2);
    seq.norm_target_sz = round([norm_target_sz_h norm_target_sz_w]);
	
    % distance (on one side) between target and bg area
	norm_pad = floor((seq.norm_bg_area - seq.norm_target_sz) / 2);
	radius = min(norm_pad);
	% norm_delta_area is the number of rectangles that are considered.
	% it is the "sampling space" and the dimension of the final merged resposne
	% it is squared to not privilege any particular direction
	seq.norm_delta_area = (2*radius+1) * [1, 1];
	% Rectangle in which the integral images are computed.
	% Grid of rectangles ( each of size norm_target_sz) has size norm_delta_area.
	seq.norm_pwp_search_area = seq.norm_target_sz + seq.norm_delta_area - 1;

end
