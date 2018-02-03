function center_likelihood = getCenterLikelihood(object_likelihood, m)
%GETCENTERLIKELIHOOD computes the sum over rectangles of size M.
% CENTER_LIKELIHOOD is the 'colour response'
    [h,w] = size(object_likelihood);
    n1 = h - m(1) + 1;
    n2 = w - m(2) + 1;

    % integral images
    SAT = integralImage(object_likelihood);
    i = 1:n1;
    j = 1:n2;
    center_likelihood = (SAT(i,j) + SAT(i+m(1), j+m(2)) - SAT(i+m(1), j) - SAT(i, j+m(2))) / prod(m);
end
