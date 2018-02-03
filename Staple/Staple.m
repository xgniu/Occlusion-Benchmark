classdef Staple < handle
    properties (SetAccess = private)
        params; % tracker parameters
        seq;    % sequence info
        scale;  %tracker.scale estimation
        tmpl;   % correlation filter template
        hist;   % color-based histogram
        state;  % state of tracking        
    end
    
    methods
        
        function tracker = Staple(im, init_pos, target_sz)
            %% Parameters setting
            tracker.params.hog_cell_size = 4;
            tracker.params.fixed_area = 150^2;           % standard area to which we resize the target
            tracker.params.n_bins = 2^5;                 % number of bins for the color histograms (bg and fg models)
            tracker.params.learning_rate_pwp = 0.04;     % bg and fg color models learning rate
            tracker.params.feature_type = 'fhog';        % 'gray';
            tracker.params.inner_padding = 0.2;          % defines inner area used to sample colors from the foreground
            tracker.params.output_sigma_factor = 0.1 ;   % standard deviation for the desired translation filter output
            tracker.params.lambda = 1e-3;                % regularization weight
            tracker.params.learning_rate_cf = 0.01;      % HOG model learning rate
            tracker.params.merge_factor = 0.3;           % fixed interpolation factor - how to linearly combine the two responses
            tracker.params.merge_method = 'const_factor';
            tracker.params.den_per_channel = false;
            
            % scale related
            tracker.params.hog_scale_cell_size = 4;
            tracker.params.learning_rate_scale = 0.025;
            tracker.params.scale_sigma_factor = 1/4;
            tracker.params.num_scales = 33;
            tracker.params.scale_model_factor = 1.0;
            tracker.params.scale_step = 1.02;
            tracker.params.scale_model_max_area = 32*16;
            
            %% Sequence info loading
            tracker.seq.curr_im = im;
            
            tracker.seq.init_pos = init_pos;
			tracker.seq.init_target_sz = target_sz;
            tracker.seq.grayscale_sequence = false;
            if(size(im,3)==1)
                tracker.seq.grayscale_sequence = true;
            end
            
            tracker.seq = initializeAllAreas(im, tracker.params, tracker.seq);
            
            %tracker.seq.videoPlayer = vision.VideoPlayer('Position', [100 100 [size(im,2), size(im,1)]+30]);
            
            %% Template initializing
            tracker.tmpl.hann_window = single(hann(tracker.seq.cf_response_size(1)) * hann(tracker.seq.cf_response_size(2))'); % 2-D
            tracker.tmpl.output_sigma = sqrt(prod(tracker.seq.norm_target_sz)) * tracker.params.output_sigma_factor / tracker.params.hog_cell_size;
            tracker.tmpl.yf = fft2( gaussianResponse(tracker.seq.cf_response_size, tracker.tmpl.output_sigma) ); % 2-D
            
            im_patch_cf = getSubwindow(im,tracker.seq.init_pos,tracker.seq.norm_bg_area, tracker.seq.bg_area);
            xt = getFeatureMap(im_patch_cf, tracker.params.feature_type, tracker.seq.cf_response_size, tracker.params.hog_cell_size);% 3-D, shape:[cf_response_sz, 28]
            xtf = fft2( bsxfun(@times, tracker.tmpl.hann_window, xt) );
            tracker.tmpl.hf_num = bsxfun(@times, conj(tracker.tmpl.yf), xtf) / prod(tracker.seq.cf_response_size); % same shape as xtf, complex
            tracker.tmpl.hf_den = (conj(xtf) .* xtf) / prod(tracker.seq.cf_response_size); %single
            
            %% Color histogram initializing
            im_patch_hist = getSubwindow(im, tracker.seq.init_pos, tracker.seq.norm_bg_area, tracker.seq.bg_area);
            tracker.hist.new_pwp_model = true;
            [tracker.hist.bg_hist, tracker.hist.fg_hist] = updateHistModel(tracker.hist.new_pwp_model, ...
                im_patch_hist, tracker.seq.bg_area, tracker.seq.fg_area, tracker.seq.init_target_sz, ...
                tracker.seq.norm_bg_area, tracker.params.n_bins, tracker.seq.grayscale_sequence);
            tracker.hist.new_pwp_model = false;
            
            %% Scale estimation
            tracker.scale.scale_factor = 1;
            tracker.scale.base_target_sz = tracker.seq.init_target_sz;
            
            tracker.scale.scale_sigma = sqrt(tracker.params.num_scales) * tracker.params.scale_sigma_factor;
            ss = (1:tracker.params.num_scales) - ceil(tracker.params.num_scales/2);
            ys = exp(-0.5 * (ss.^2) / tracker.scale.scale_sigma^2);
            tracker.scale.ysf = single(fft(ys));
            
            tracker.scale.scale_window = single(hann(tracker.params.num_scales));
                       
            ss = 1:tracker.params.num_scales;
            tracker.scale.scale_factors = tracker.params.scale_step.^(ceil(tracker.params.num_scales/2) - ss);
            
            if tracker.params.scale_model_factor^2 * prod(tracker.seq.norm_target_sz) > tracker.params.scale_model_max_area
                tracker.params.scale_model_factor = sqrt(tracker.params.scale_model_max_area/prod(tracker.seq.norm_target_sz));
            end
            
            tracker.scale.scale_model_sz = floor(tracker.seq.norm_target_sz * tracker.params.scale_model_factor);
            
            tracker.scale.min_scale_factor = tracker.params.scale_step ^ ceil(log(max(5 ./ tracker.seq.bg_area)) ...
                / log(tracker.params.scale_step));
            tracker.scale.max_scale_factor = tracker.params.scale_step ^ floor(log(min([size(im,1) size(im,2)] ...
                ./ tracker.seq.init_target_sz)) / log(tracker.params.scale_step));
            
            im_patch_scale = getScaleSubwindow(im, tracker.seq.init_pos, tracker.scale.base_target_sz, ...
                tracker.scale.scale_factor*tracker.scale.scale_factors, tracker.scale.scale_window, ...
                tracker.scale.scale_model_sz, tracker.params.hog_scale_cell_size);
            xsf = fft(im_patch_scale,[],2);
            tracker.scale.sf_num = bsxfun(@times, tracker.scale.ysf, conj(xsf));
            tracker.scale.sf_den = sum(xsf .* conj(xsf), 1);
            
            %% Tracking state
            tracker.state.center = tracker.seq.init_pos;
            tracker.state.target_sz = tracker.seq.init_target_sz;
            tracker.state.rect_position = [tracker.seq.init_pos([2,1]) - tracker.seq.init_target_sz([2,1])/2, tracker.seq.init_target_sz([2,1])];
            tracker.state.rect_position_padded = [tracker.seq.init_pos([2,1]) - tracker.seq.bg_area([2,1])/2, tracker.seq.bg_area([2,1])];
            
            %im_show = insertShape(im, 'Rectangle', tracker.state.rect_position, 'LineWidth', 4, 'Color', 'black');
            %im_show = insertShape(im_show, 'Rectangle', tracker.state.rect_position_padded, 'LineWidth', 4, 'Color', 'yellow');
            %step(tracker.seq.videoPlayer, im_show);
        end
        
        function [track_output, response] = track(self,im)
            %% Tracking based on Template 
            % extract patch of size bg_area and resize to norm_bg_area
            self.seq.curr_im = im;
            im_patch_cf = getSubwindow(im, self.state.center, self.seq.norm_bg_area, self.seq.bg_area);
            xt = getFeatureMap(im_patch_cf, self.params.feature_type, self.seq.cf_response_size, self.params.hog_cell_size);
            xtf = fft2( bsxfun(@times, self.tmpl.hann_window, xt) );
            
            if self.params.den_per_channel
                hf = self.tmpl.hf_num ./ (self.tmpl.hf_den + self.params.lambda);
            else
                hf = bsxfun(@rdivide, self.tmpl.hf_num, sum(self.tmpl.hf_den, 3)+self.params.lambda);
            end
            response_cf = self.ensure_real( ifft2(sum(conj(hf) .* xtf, 3)) );
            
            % Crop square search region (in feature pixels).
            response_cf = cropFilterResponse(response_cf, self.floor_odd(self.seq.norm_delta_area / self.params.hog_cell_size));
            if self.params.hog_cell_size > 1
                % Scale up to match center likelihood resolution.
                response_cf = mexResize(response_cf, self.seq.norm_delta_area,'auto');
            end
            
            %% Tracking based on Histogram
            pwp_search_area = round(self.seq.norm_pwp_search_area / self.seq.area_resize_factor);
            % extract patch of size pwp_search_area and resize to norm_pwp_search_area
            im_patch_pwp = getSubwindow(im, self.state.center, self.seq.norm_pwp_search_area, pwp_search_area);
            likelihood_map = getColourMap(im_patch_pwp, self.hist.bg_hist, self.hist.fg_hist, self.params.n_bins, self.seq.grayscale_sequence);
            likelihood_map(isnan(likelihood_map)) = 0;
            
            % each pixel of response_pwp loosely represents the likelihood that
            % the target (of size norm_target_sz) is centred on it
            response_pwp = getCenterLikelihood(likelihood_map, self.seq.norm_target_sz);
            
            %% Merge response and find the maximun response. Position estimation
            response = mergeResponses(response_cf, response_pwp, self.params.merge_factor, self.params.merge_method);
            [row, col] = find(response == max(response(:)), 1);
            center = (1+self.seq.norm_delta_area) / 2;
            self.state.center = self.state.center + ([row, col] - center) / self.seq.area_resize_factor;
            
            %% Scale estimation
            im_patch_scale = getScaleSubwindow(im, self.state.center, self.scale.base_target_sz, ...
                self.scale.scale_factor * self.scale.scale_factors, self.scale.scale_window,...
                self.scale.scale_model_sz, self.params.hog_scale_cell_size);
            xsf = fft(im_patch_scale,[],2);
            scale_response = real(ifft(sum(self.scale.sf_num .* xsf, 1) ./ (self.scale.sf_den + self.params.lambda) ));
            recovered_scale = ind2sub(size(scale_response),find(scale_response == max(scale_response(:)), 1));
            self.scale.scale_factor = self.scale.scale_factor * self.scale.scale_factors(recovered_scale);
            
            if self.scale.scale_factor < self.scale.min_scale_factor
                self.scale.scale_factor = self.scale.min_scale_factor;
            elseif self.scale.scale_factor > self.scale.max_scale_factor
                self.scale.scale_factor = self.scale.max_scale_factor;
            end
            
            self.state.target_sz = round(self.scale.base_target_sz * self.scale.scale_factor);
            
            track_output = [self.state.center([2,1]) - self.state.target_sz([2,1])/2, self.state.target_sz([2,1])];
            
            self.state.rect_position = track_output;
            self.state.rect_position_padded = [self.state.center([2,1]) - self.seq.bg_area([2,1])/2, self.seq.bg_area([2,1])];
            
            %im_show = insertShape(im, 'Rectangle', self.state.rect_position, 'LineWidth', 4, 'Color', 'black');
            %im_show = insertShape(im_show, 'Rectangle', self.state.rect_position_padded, 'LineWidth', 4, 'Color', 'yellow');
            %step(self.seq.videoPlayer, im_show);
        end
        
        function update(self)
            %% propertity 'seq' update
            im = self.seq.curr_im;
            avg_dim = sum(self.state.target_sz)/2;
            bg_area = round(self.state.target_sz + avg_dim);
            if(bg_area(2)>size(im,2)),  bg_area(2)=size(im,2)-1;    end
            if(bg_area(1)>size(im,1)),  bg_area(1)=size(im,1)-1;    end           
            self.seq.bg_area = bg_area - mod(bg_area - self.state.target_sz, 2);
            
            fg_area = round(self.state.target_sz - avg_dim * self.params.inner_padding);
            self.seq.fg_area = fg_area + mod(self.seq.bg_area - fg_area, 2);
            
            self.seq.area_resize_factor = sqrt(self.params.fixed_area/prod(self.seq.bg_area));
            
            %% template update
            im_patch_cf = getSubwindow(im, self.state.center, self.seq.norm_bg_area, self.seq.bg_area);
            xt = getFeatureMap(im_patch_cf, self.params.feature_type, self.seq.cf_response_size, self.params.hog_cell_size);% 3-D, shape:[cf_response_sz, 28]
            xtf = fft2( bsxfun(@times, self.tmpl.hann_window, xt) );
            
            new_hf_num = bsxfun(@times, conj(self.tmpl.yf), xtf) / prod(self.seq.cf_response_size); % same shape as xtf, complex
            new_hf_den = (conj(xtf) .* xtf) / prod(self.seq.cf_response_size); %single

            self.tmpl.hf_den = (1 - self.params.learning_rate_cf) * self.tmpl.hf_den + self.params.learning_rate_cf * new_hf_den;
            self.tmpl.hf_num = (1 - self.params.learning_rate_cf) * self.tmpl.hf_num + self.params.learning_rate_cf * new_hf_num;
            
            %% hist update
            [self.hist.bg_hist, self.hist.fg_hist] = updateHistModel(self.hist.new_pwp_model, im_patch_cf,...
                self.seq.bg_area, self.seq.fg_area, self.state.target_sz, self.seq.norm_bg_area, ...
                self.params.n_bins, self.seq.grayscale_sequence, self.hist.bg_hist, self.hist.fg_hist, self.params.learning_rate_pwp);
            
            %% scale update
            im_patch_scale = getScaleSubwindow(im, self.state.center, self.scale.base_target_sz, ...
                self.scale.scale_factor * self.scale.scale_factors, self.scale.scale_window,...
                self.scale.scale_model_sz, self.params.hog_scale_cell_size);
            xsf = fft(im_patch_scale,[],2);
            new_sf_num = bsxfun(@times, self.scale.ysf, conj(xsf));
            new_sf_den = sum(xsf .* conj(xsf), 1);
            self.scale.sf_den = (1 - self.params.learning_rate_scale) * self.scale.sf_den + self.params.learning_rate_scale * new_sf_den;
            self.scale.sf_num = (1 - self.params.learning_rate_scale) * self.scale.sf_num + self.params.learning_rate_scale * new_sf_num;            
            
        end
        
        function y = floor_odd(self,x)
            y = 2*floor((x-1) / 2) + 1;
        end
        
        function y = ensure_real(self,x)
            assert(norm(imag(x(:))) <= 1e-5 * norm(real(x(:))));
            y = real(x);
        end
    end
end