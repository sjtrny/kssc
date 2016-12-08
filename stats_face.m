paths = ['common:'];
addpath(paths);

run('vlfeat/toolbox/vl_setup');

if exist('matlabpool','file') && matlabpool('size') == 0
    matlabpool open
end

rng(1);

n_pic = 64;
% n_pic = 15;

class_inds = [1:13, 15:39];

n_class = size(class_inds, 2);

n_samples = n_pic * n_class;

X = zeros(32256, n_samples);

% for k = 1 : n_class
%     
%     k_list = dir(['data/CroppedYale/yaleB' num2str(class_inds(k), '%02.0f')  '/*P00A*.pgm*']);
%     
%     for j = 1 : n_pic
%         
%         try
%         X(:, k*n_pic - n_pic + j) = reshape(double(imread(['data/CroppedYale/yaleB' num2str(class_inds(k), '%02.0f') '/' k_list(j).name]))/255, 32256, 1);
%         catch
%             disp('nothing');
%         end
%         
%     end
%     
% end

tic;
for k = 1 : n_class
    
    k_list = dir(['data/CroppedYale/yaleB' num2str(class_inds(k), '%02.0f')  '/*P00A*.pgm*']);
    
    ref_image = double(imread(['data/CroppedYale/yaleB' num2str(class_inds(k), '%02.0f') '/' k_list(1).name]))/255;
    
    for j = 1 : n_pic
        
        try
        X(:, k*n_pic - n_pic + j) = reshape( imhistmatch( double(imread(['data/CroppedYale/yaleB' num2str(class_inds(k), '%02.0f') '/' k_list(j).name]))/255, ref_image), 32256, 1);
        catch
            disp('nothing');
        end
        
    end
    
end
toc;

tic;
true_pos = get_block_diag(n_samples, n_class);
false_pos = 1 - true_pos;

max_nn = 30;

true_pos_rate = zeros(max_nn - 1,1);
false_pos_rate = zeros(max_nn - 1,1);

vlfeat_tree = vl_kdtreebuild(X);

parfor k = 1 : max_nn;
    
    full_rows_mat = double(vl_kdtreequery(vlfeat_tree,X,X,'NUMNEIGHBORS',k+1,'MAXNUMCOMPARISONS',10*(k+1/10)));
    
    rows_mat = reshape(full_rows_mat(2:end, :), k*n_samples, 1);
    cols_mat = reshape(repmat(1:n_samples, k, 1), k*n_samples, 1);
    
    index_mat = sparse(rows_mat, cols_mat, ones(k*n_samples, 1), n_samples, n_samples);
    
    true_pos_rate(k) = nnz(index_mat .* true_pos);
    false_pos_rate(k) = nnz(index_mat .* false_pos);
    
end
toc;

tic;
singular_values = zeros(n_pic, n_class);
for k = 1 : n_class
    
    [U, S, V] = svds(X(:, 1 + ((k-1)*n_pic) : k*n_pic ), n_pic);
    
    singular_values(:, k) = diag(S);
    
end
toc;

plot(singular_values(1:30,1:5))

set(gca, 'fontsize', 14);
xlabel('Singular Value Indices', 'FontSize', 18);
ylabel('Singular Values', 'FontSize', 18);

print(gcf, '-depsc2', 'singular_values_yale.eps');

close all

normalized_true_pos = (true_pos_rate ./ ((1:(max_nn))'*n_samples)) * 100;
normalized_false_pos = (false_pos_rate ./ ((1:(max_nn))'*n_samples)) * 100;



figure, h1 = plot(normalized_true_pos, 'b');
hold on
h2 = plot(normalized_false_pos, 'r');

set(gca,'YTickLabel', {'0', '20', '40', '60', '80', '100', ''})



h3 = plot(enough_true_pos, 'g');

set(gca, 'fontsize', 14);


xlabel('k - Number of Nearest Neighbours', 'FontSize', 18);
ylabel('Percent %', 'FontSize', 18);

xlim([1, 30]);
ylim([-5, 120]);

% left_axes = gca;
% 
% axesPosition = get(gca,'Position');          %# Get the current axes position
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

% enough_true_pos = ((true_pos_rate ./ (n_samples)) >= 9)*100;

ylim([-5, 120]);

legend([h1, h2, h3], 'True Positives', 'False Positives', 'Sufficient True Positives', 'Location', 'NorthWest');

print(gcf, '-depsc2', 'nearest_k_yale.eps');


rmpath(paths);