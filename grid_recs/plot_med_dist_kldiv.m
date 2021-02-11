[gbins, gedges] = discretize(grid_distances, 10);
grid_med_kldivs = zeros(1, 10);
for kB = 1:10
grid_med_kldivs(kB) = median(grid_kl_divs(gbins == kB));
end
gcenters = mean([gedges(1:end-1); gedges(2:end)]);

pedges = -70:140:1330;
[pbins] = discretize(probe_distances, pedges);
pcenters = mean([pedges(1:end-1); pedges(2:end)]);
probe_med_kldivs = zeros(1, 10);

for kB = 1:10
probe_med_kldivs(kB) = median(probe_kl_divs(pbins == kB));
end

hold on;
plot(pcenters, probe_med_kldivs, 'b*-', 'LineWidth', 2, 'DisplayName', 'Median');
plot(gcenters, grid_med_kldivs, 'b*-', 'LineWidth', 2, 'DisplayName', 'Median');
