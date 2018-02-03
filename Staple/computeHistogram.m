function histogram = computeHistogram(patch, mask, n_bins, grayscale_sequence)
%COMPUTEHISTOGRAM creates a colour (or grayscale) histogram of an image patch
% MASK has the same size as the image patch and selects what should
% be used when computing the histogram (i.e. out-of-frame regions are ignored)

	[h, w, d] = size(patch);

	assert(all([h w]==size(mask)) == 1, 'mask and image are not the same size');

	bin_width = 256/n_bins; % n_bins = 32. bin_width = 8

	% convert image to 1d array with same n channels of img patch
	patch_array = reshape(double(patch), w*h, d);
	% compute to which bin each pixel (for all 3 channels) belongs to
	bin_indices = floor(patch_array/bin_width) + 1; % from 1 to 33

	if grayscale_sequence
        % A = accumarray(subs,val,sz) : A(i)=sum( val(subs(:)==i)) ). size(A): sz.
		histogram = accumarray(bin_indices, mask(:), [n_bins 1]) / sum(mask(:));
	else
		% the histogram is a cube of side n_bins
		histogram = accumarray(bin_indices, mask(:), [n_bins n_bins n_bins]) / sum(mask(:)); %size:[32,32,32]
	end

end
% classify the grayscale values into groups from 1 to 33.
% count the frequency of each group: #(pixels of the group) / #(pixels of bg or fg).