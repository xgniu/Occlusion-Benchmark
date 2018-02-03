classdef DSST < handle
    % The following properties can be set only by class methods
    properties (SetAccess = private)
        padding = 1.5;      % extra area surrounding the target
        lambda = 1e-3;  	% regularization weight (denoted "lambda" in the paper)
        interp_factor = 0.02; % tracking model learning rate (denoted "eta" in the paper)
        learning_rate = 0.02; % for scale filter
        
        output_sigma_factor = 0.1; 	% standard deviation for the desired translation filter output
        output_sigma;
        scale_sigma_factor = 0.25;        % standard deviation for the desired scale filter output
        scale_sigma;
        kernel_sigma=0.5;
        
        hog_cell_sz=4;
        hog_orientations=9;
        
        nScales = 33;         % number of scale levels (denoted "S" in the paper)
        scale_step = 1.02;               % Scale increment factor (denoted "a" in the paper)
        scale_model_max_area = 1024;      % the maximum size of scale examples
        scaleFactors;
        currentScaleFactor;
        min_scale_factor;
        max_scale_factor;
        % desired scale filter output (gaussian shaped), bandwidth proportional to number of scales
        scale_model_sz;
        
        init_target_sz;
        base_target_sz;
        center;
        target_sz;
        window_sz;
        box;
        
        cos_window;
        scale_window;
        
        model_alphaf;
        model_xf;
        sf_den;
        sf_num;
        y;
        yf;
        ysf;
        current_img;

    end
    
    methods
        function tracker = DSST(im, center, target_sz)
            
            tracker.window_sz = floor(target_sz * (1 + tracker.padding));
            tracker.center = center;
            tracker.target_sz = target_sz;
            tracker.init_target_sz = target_sz;
            % target size at scale = 1
            tracker.base_target_sz = target_sz;
                       
            tracker.box = [ center([2 1]) - target_sz([2 1])/2, target_sz([2 1]) ];%(x,y,width,height)
            
            % desired scale filter output (gaussian shaped), bandwidth proportional to number of scales
            tracker.scale_sigma = sqrt(tracker.nScales) * tracker.scale_sigma_factor;%nScales=33. scale_sigma_factor=1/4.
            ss = (1:tracker.nScales) - ceil(tracker.nScales/2);%(1:33)-17=-16:16
            ys = exp(-0.5 * (ss.^2) / tracker.scale_sigma^2);%1-D
            tracker.ysf = single(fft(ys));
            
            if mod(tracker.nScales,2) == 0
                tracker.scale_window = single(hann(tracker.nScales+1));
                tracker.scale_window = tracker.scale_window(2:end);
            else
                tracker.scale_window = single(hann(tracker.nScales));
            end
            
            % scale factors
            ss = 1:tracker.nScales;%1:33
            tracker.scaleFactors = tracker.scale_step.^( ceil(tracker.nScales/2) - ss );%scale_step=1.02. 17-1:33=16:-16
            
            % compute the resize dimensions used for feature extraction in the scale estimation
            scale_model_factor = 1;
            if prod(tracker.init_target_sz) > tracker.scale_model_max_area
                scale_model_factor = sqrt(tracker.scale_model_max_area/prod(tracker.init_target_sz));
            end
            tracker.scale_model_sz = floor(tracker.init_target_sz * scale_model_factor);%constant
            
            tracker.currentScaleFactor = 1;
            tracker.min_scale_factor = tracker.scale_step ^ ceil(log(max(5 ./ tracker.window_sz)) / log(tracker.scale_step));
            tracker.max_scale_factor = tracker.scale_step ^ floor(log(min([size(im,1) size(im,2)] ./ tracker.base_target_sz)) / log(tracker.scale_step));
            
            %create regression labels, gaussian shaped, with a bandwidth proportional to target size
            tracker.output_sigma = sqrt(prod(tracker.target_sz)) * tracker.output_sigma_factor;
            tracker.y = gaussian_shaped_labels(tracker.output_sigma, floor(tracker.window_sz) );
            tracker.yf = single( fft2( tracker.y ) );%constant
            
            %store pre-computed cosine window
            tracker.cos_window = hann(size(tracker.yf,1)) * hann(size(tracker.yf,2))';%constant
            
            patch = ICF_get_subwindow(im, center, tracker.window_sz, tracker.currentScaleFactor);
            xf = fft2(ICF_get_features(patch, tracker.cos_window));
            tracker.model_xf = xf;
            
            %calculate response of the classifier at all shifts
            kf = gaussian_correlation(xf, xf, tracker.kernel_sigma);

            alphaf = tracker.yf ./ (kf + tracker.lambda);   %equation for fast training
            tracker.model_alphaf = alphaf;
            
            % extract the training sample feature map for the scale filter
            xs = get_scale_sample(im, center, tracker.base_target_sz,...
                tracker.currentScaleFactor * tracker.scaleFactors, tracker.scale_window, tracker.scale_model_sz);
            % calculate the scale filter update
            xsf = fft(xs,[],2);
            new_sf_num = bsxfun(@times, tracker.ysf, conj(xsf));
            new_sf_den = sum(xsf .* conj(xsf), 1);
            tracker.sf_num = new_sf_num;
            tracker.sf_den = new_sf_den;
            
        end
        
        function [box, response] = track(self, im)
            self.current_img = im;
            patch = ICF_get_subwindow(im, self.center, self.window_sz, self.currentScaleFactor);
            zf = fft2(ICF_get_features(patch, self.cos_window));

            kzf = gaussian_correlation(zf, self.model_xf, self.kernel_sigma);
            
            response =sum( real( ifft2(self.model_alphaf .* kzf) ), 3 );%equation for fast detection
            [vert_delta, horiz_delta] = find(response == max(response(:)), 1);
            if vert_delta > size(zf,1) / 2  %wrap around to negative half-space of vertical axis
                vert_delta = vert_delta - size(zf,1);
            end
            if horiz_delta > size(zf,2) / 2  %same for horizontal axis
                horiz_delta = horiz_delta - size(zf,2);
            end
            center1 = self.center + round([vert_delta - 1, horiz_delta - 1] * self.currentScaleFactor);
            if all(center1>1) && all( ([size(im,1) size(im,2)]-center1)>1 )
                self.center=center1;
            end
            
            xs = get_scale_sample(im, self.center, self.base_target_sz,...
                self.currentScaleFactor * self.scaleFactors, self.scale_window, self.scale_model_sz);
            xsf = fft(xs,[],2);
            scale_response = real(ifft(sum(self.sf_num .* xsf, 1) ./ (self.sf_den + self.lambda)));
            
            % find the maximum scale response
            recovered_scale = find(scale_response == max(scale_response(:)), 1);
            % update the scale
            self.currentScaleFactor = self.currentScaleFactor * self.scaleFactors(recovered_scale);
            if self.currentScaleFactor < self.min_scale_factor
                self.currentScaleFactor = self.min_scale_factor;
            elseif self.currentScaleFactor > self.max_scale_factor
                self.currentScaleFactor = self.max_scale_factor;
            end
            
            self.target_sz = floor(self.base_target_sz * self.currentScaleFactor);
            
            position = self.center([2 1]) - self.target_sz([2 1])/2;
            box = [position, self.target_sz([2 1])];
            self.box=box;
        end
        
        function update(self)
            cur_interp_factor = self.interp_factor;
            cur_learning_rate = self.learning_rate;
            
            im = self.current_img;
            %obtain a subwindow for training at newly estimated target position
            patch = ICF_get_subwindow(im, self.center, self.window_sz, self.currentScaleFactor);            
            xf = fft2(ICF_get_features(patch, self.cos_window));
            self.model_xf = (1 - cur_interp_factor) * self.model_xf + cur_interp_factor * xf;
            
            kf = gaussian_correlation(xf, xf, self.kernel_sigma);

            alphaf = self.yf ./ (kf + self.lambda);   %equation for fast training
            self.model_alphaf = (1 - cur_interp_factor) * self.model_alphaf + cur_interp_factor * alphaf;
            
            % extract the training sample feature map for the scale filter
            xs = get_scale_sample(im, self.center, self.base_target_sz, ...
                self.currentScaleFactor * self.scaleFactors, self.scale_window, self.scale_model_sz);
            % calculate the scale filter update
            xsf = fft(xs,[],2);
            new_sf_num = bsxfun(@times, self.ysf, conj(xsf));
            new_sf_den = sum(xsf .* conj(xsf), 1);
            self.sf_den = (1 - cur_learning_rate) * self.sf_den + cur_learning_rate * new_sf_den;
            self.sf_num = (1 - cur_learning_rate) * self.sf_num + cur_learning_rate * new_sf_num;
        end
    end % methods
end % classdef