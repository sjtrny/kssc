function [ labels ] = get_sssc_clusters( X, n_training, lambda, n_clusters )
%GET_SSSC_CLUSTERS Summary of this function goes here
%   Detailed explanation goes here

[~, n] = size(X);

tol = 1*10^-6;
par.maxIteration = 5000;
par.isNonnegative = false;

training_inds = randsample(n, n_training);

out_inds = 1:n;
out_inds(training_inds) = [];

Tr_X = X(:, training_inds);
Out_X = X(:, out_inds);

training_labels = InSample(Tr_X, lambda, tol, par, n_clusters);

out_labels = OutSample(Tr_X, Out_X, training_labels);

labels = zeros(n, 1);
labels(training_inds) = training_labels;
labels(out_inds) = out_labels;

end

