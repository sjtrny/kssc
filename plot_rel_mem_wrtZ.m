resolution_steps = 10;

max_N = 100000;

k_min = 1;
k_max = max_N;

D = 500;

kssc_exact = zeros(resolution_steps + 1, 1);
kssc_relaxed = zeros(resolution_steps + 1, 1);

k_list = zeros(resolution_steps + 1, 1);

for k = 1 : resolution_steps + 1;
    
    current_k = ((k_max - k_min) / resolution_steps) * (k-1);
    
    k_list(k) = current_k;
    
    kssc_exact(k) = 2*current_k*max_N + D*max_N;
    kssc_relaxed(k) = 4*current_k*max_N;
    
end

ssc_exact = repmat(2 * max_N^2 + D*max_N, resolution_steps + 1, 1);
ssc_exact = (ssc_exact * 64) * 1.25e-10;

ssc_relaxed = repmat(4*max_N^2, resolution_steps + 1, 1);
ssc_relaxed = (ssc_relaxed * 64) * 1.25e-10;

kssc_exact = (kssc_exact * 64) * 1.25e-10;
kssc_relaxed = (kssc_relaxed * 64) * 1.25e-10;

h1 = plot(ssc_exact, '-*b');
hold on
h2 = plot(ssc_relaxed, '-ob');

h3 = plot(kssc_exact, '-*r');
h4 = plot(kssc_relaxed, '-or');

set(gca, 'fontsize', 14);

xlabel('k as percentage of N', 'FontSize', 18);
xlim([1, 11]);
set(gca, 'XTick', 1:2:(resolution_steps + 1));
set(gca, 'XTickLabel', round(100*k_list(1:2:length(k_list))/max_N));

ylabel('RAM Required wrt Z - Gigabytes', 'FontSize', 18);
ylim([-10, 350]);

legend([h1, h2, h3, h4], 'SSC Exact', 'SSC Relaxed', 'kSSC Exact', 'kSSC Relaxed', 'Location', 'SouthEast');

print(gcf, '-depsc2', 'relmem.eps');
