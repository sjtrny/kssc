function [ rows, cols ] = process_neighbours( rows_mat )
% Removes diagonal entries and spits out clean vectors of row/col indices

[m, n] = size(rows_mat);

cols = zeros(m * n, 1);

rows = zeros(m * n, 1);

count = 1;

for k = 1 : n
    sub_rows = rows_mat(rows_mat(:, k) ~= k, k);
    if (length(sub_rows) >= m)
        sub_rows = sub_rows(1:m-1);
    end
    rows(count : count + length(sub_rows) - 1) = sub_rows;
    cols(count : count + length(sub_rows) - 1) = k;
    count = count + length(sub_rows);
end

rows = rows(1:count - 1);
cols = cols(1:count - 1);

end

