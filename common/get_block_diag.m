function [ X ] = get_block_diag(dim, n_blocks)
%GET_BLOCK_DIAG Summary of this function goes here
%   Detailed explanation goes here

if ~exist('dim', 'var')
    error('Matrix size not provided.');
end

if ~exist('n_blocks', 'var')
    error('Number of blocks not provided.');
end

if mod(dim, n_blocks)
    error('n_blocks is not a factor of dim');
end

block_size = dim / n_blocks;

X = zeros(dim);

for k = 1 : n_blocks
    X(1 + (k-1)*block_size: k*block_size, 1 + (k-1)*block_size: k*block_size) = 1;
end

end

