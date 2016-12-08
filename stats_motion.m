paths = ['common:'];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end

rng(1);

% SVD each motion
% For each capture do nearest neighbour false/true pos

file_list = dir('data/Hopkins155/*.mat');

n_captures = length(file_list);

n_singular_values = 30;

singular_values = zeros(30, 1);

max_nn = 30;

normalized_true_pos = zeros(max_nn, n_captures);
normalized_false_pos = zeros(max_nn, n_captures);

enough_true_pos = zeros(max_nn, n_captures);

for k = 1 : n_captures
    
    load(['data/Hopkins155/' file_list(k).name]);
    
    n_class = length(unique(s));
    n_samples = points;

    X = zeros(2*frames, points);
    for f = 1 : frames
        X(f*2 - 1 : f*2, :) = x(1:2, :, f);
    end

    [a, b] = sort(s);

    X = X(:,b);
    
    for j = 1 : n_class
        
        [U, S, V] = svds(X(:, a == j), n_singular_values);

        max_singular_value = min(n_singular_values, max(size(S)));
        
        singular_values(1:max_singular_value, k) = diag(S);
        
    end   

    true_pos = zeros(n_samples);
    position = 1;
    for f = 1 : n_class
        n_in_class = length(find(s == f));
        true_pos(position:(position-1)+n_in_class, position:(position-1)+n_in_class) = 1;
        position = position + n_in_class;
    end
    
    false_pos = 1 - true_pos;


    true_pos_rate = zeros(max_nn,1);
    false_pos_rate = zeros(max_nn ,1);

    vlfeat_tree = vl_kdtreebuild(X);
    
    for j = 1 : max_nn;
    
        full_rows_mat = double(vl_kdtreequery(vlfeat_tree,X,X,'NUMNEIGHBORS',j+1,'MAXNUMCOMPARISONS',10*(j+1/10)));

        rows_mat = reshape(full_rows_mat(2:end, :), j*n_samples, 1);
        cols_mat = reshape(repmat(1:n_samples, j, 1), j*n_samples, 1);
        
        index_mat = sparse(rows_mat, cols_mat, ones(j*n_samples, 1), n_samples, n_samples);
        
        true_pos_rate(j) = nnz(index_mat .* true_pos);
        false_pos_rate(j) = nnz(index_mat .* false_pos);
        
    
    end

    normalized_true_pos(:,k) = (true_pos_rate ./ ((1:(max_nn))'*n_samples) ) * 100;
    normalized_false_pos(:,k) = (false_pos_rate ./ ((1:(max_nn))'*n_samples) ) * 100;
    
    enough_true_pos(:,k) = (true_pos_rate ./ (n_samples)) >= 4;
    
end

enough_true_pos = enough_true_pos * 100;


plot(singular_values(1:30,1:5))

set(gca, 'fontsize', 14);
xlabel('Singular Value Indices', 'FontSize', 18);
ylabel('Singular Values', 'FontSize', 18);

print(gcf, '-depsc2', 'singular_values_hopkins.eps');

close all



figure, h1 = plot(mean(normalized_true_pos,2), 'b');
hold on
h2 = plot(mean(normalized_false_pos,2), 'r');

h3 = plot(median(enough_true_pos,2), 'g');

set(gca, 'fontsize', 14);


xlabel('k - Number of Nearest Neighbours', 'FontSize', 18);
ylabel('Percent %', 'FontSize', 18);

xlim([1, 30]);
ylim([-5, 140]);


set(gca,'YTickLabel', {'0', '20', '40', '60', '80', '100', '', '', ''})

% left_axes = gca;
% 
axesPosition = get(gca,'Position');          %# Get the current axes position
hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
                'Color','none',...           %#   ... with no background color
                'YLim',[0 1],...            %#   ... and a different scale
                'YAxisLocation','right',...  %#   ... located on the right
                'XTick',[],...               %#   ... with no x tick marks
                'Box','off');                %#   ... and no surrounding box
set(hNewAxes, 'fontsize', 14);
ylabel('Sufficient True Positives', 'FontSize', 18);
set(hNewAxes,'YTick',[0, 100])
set(hNewAxes,'YTickLabel', {'False', 'True'})

ylim([-5, 140]);

legend([h1, h2, h3], 'True Positives', 'False Positives', 'Sufficient True Positives', 'Location', 'NorthWest');

print(gcf, '-depsc2', 'nearest_k_hopkins.eps');




rmpath(paths);