resolution_steps = 10;

max_N = 100000;
D = 500;
n_k_1 = 100;
n_k_2 = 1000;

ssc_exact = zeros(resolution_steps + 1, 1);
ssc_relaxed = zeros(resolution_steps + 1, 1);

kssc_exact_1 = zeros(resolution_steps + 1, 1);
kssc_relaxed_1 = zeros(resolution_steps + 1, 1);

kssc_exact_2 = zeros(resolution_steps + 1, 1);
kssc_relaxed_2 = zeros(resolution_steps + 1, 1);

n_list = zeros(resolution_steps + 1, 1);

for k = 1 : resolution_steps + 1;
    
    current_N = (max_N / resolution_steps) * (k-1);
    
    n_list(k) = current_N;
    
    ssc_exact(k) = 3*(current_N^2) + 4*(current_N^2)*D + 4*D*current_N + 4*current_N^2;
    ssc_relaxed(k) = 2*(current_N^2) + 4*(current_N^2)*D + D*current_N + 4*current_N^2;
    
    kssc_exact_1(k) = 3*n_k_1*current_N + 4*n_k_1*current_N*D + 4*D*current_N + 4*n_k_1*current_N;
    kssc_relaxed_1(k) = 2*n_k_1*current_N + 4*n_k_1*current_N*D + D*current_N + 4*n_k_1*current_N;
    
    kssc_exact_2(k) = 3*n_k_2*current_N + 4*n_k_2*current_N*D + 4*D*current_N + 4*n_k_2*current_N;
    kssc_relaxed_2(k) = 2*n_k_2*current_N + 4*n_k_2*current_N*D + D*current_N + 4*n_k_2*current_N;
end


h1 = plot(ssc_exact, '-*b');
hold on
h2 = plot(ssc_relaxed, '-ob');

h3 = plot(kssc_exact_1, '-*r');
h4 = plot(kssc_relaxed_1, '-or');

h5 = plot(kssc_exact_2, '-*k');
h6 = plot(kssc_relaxed_2, '-ok');

set(gca, 'fontsize', 14);

xlabel('N - Number of data points', 'FontSize', 18);
xlim([1, 11]);
set(gca, 'XTick', 1:2:(resolution_steps + 1));
set(gca, 'XTickLabel', n_list(1:2:length(n_list)));

ylabel('FLOP Count', 'FontSize', 18);
init_ylabel = ylim;
init_ylabel(1) = -1000000000000;
ylim(init_ylabel);

legend([h1, h2, h3, h4, h5, h6], 'SSC Exact', 'SSC Relaxed', ['kSSC Exact (k = ' num2str(n_k_1) ')'],['kSSC Relaxed (k = ' num2str(n_k_1) ')'], ['kSSC Exact (k = ' num2str(n_k_2) ')'],['kSSC Relaxed (k = ' num2str(n_k_2) ')'], 'Location', 'NorthWest');

print(gcf, '-depsc2', 'flops.eps');

close all