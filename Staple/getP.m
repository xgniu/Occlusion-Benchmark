function P = getP(histogram, h, w, bin_indices, grayscale_sequence)
    %GETP computes the scores given the histogram
    % query the hist for the probability of each pixel
    if grayscale_sequence
        hist_indices = bin_indices;
    else
        hist_indices = sub2ind(size(histogram), bin_indices(:,1), bin_indices(:,2), bin_indices(:,3));
    end
    % hist_indices : to which bin each pixel belongs
    
    % shape it as a matrix
    P = reshape(histogram(hist_indices), h, w);
    % P : the frequency of the bin to which each pixel belongs
end 