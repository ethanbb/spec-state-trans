% Purpose: produce figures that demonstrate some properties of mutual information.

n_classes = 4;
n_time = 200;
trans_prob = 0.2;
rng(2109876308);

% Generate first channel states by a Poisson process
chan1 = ones(1, n_time);
for kT = 2:n_time
    if rand < trans_prob
        chan1(kT:end) = randsample(setdiff(1:n_classes, chan1(kT)), 1);
    end
end

stay_prob = 0.1;
look_std = 1;

%% Case 1: channels are synchronized, with noise
chan2_sync = chan1(1) * ones(size(chan1));
for kT = 1:n_time
    if rand > stay_prob
        look_offset = round(look_std * randn);
        chan2_sync(kT:end) = chan1(max(1, min(n_time, kT + look_offset)));
    end
end

plot_mutinfo_demo(chan1, chan2_sync, 'Synchronized');

%% Case 2: synchronized but class class assignments are different

chan2_iso = chan2_sync;
chan2_iso = mod(chan2_iso, n_classes) + 1;

plot_mutinfo_demo(chan1, chan2_iso, 'Isomorphic');

%% Case 3: channels are synchronized but channel 2 has a smaller repertoire

chan2_compressed = chan2_sync;
chan2_compressed(chan2_compressed == 3) = 1;
chan2_compressed(chan2_compressed == 4) = 2;

plot_mutinfo_demo(chan1, chan2_compressed, 'Compressed');

%% Case 4: channels are independent
chan2_indep = ones(1, n_time);
for kT = 2:n_time
    if rand < trans_prob
        chan2_indep(kT:end) = randsample(setdiff(1:n_classes, chan2_indep(kT)), 1);
    end
end

plot_mutinfo_demo(chan1, chan2_indep, 'Independent');

%%
function plot_mutinfo_demo(chan1, chan2, name)
figure('Position', [100, 100, 660, 170]);
ax = axes;
class_plot(1:length(chan1), [0.5, 1.5], chan1);
hold on;
class_plot(1:length(chan2), [1.5, 2.5], chan2);

ax.YDir = 'reverse';
yticks([1, 2]);
yticklabels({'\fontsize{16}\bf{A}', '\fontsize{16}\bf{B}'});

mi = class_mut_info([chan1(:), chan2(:)]);
title(sprintf('%s: mutual information = %.2f bits', name, mi(2,1)), 'FontSize', 18);

legend('off');

end