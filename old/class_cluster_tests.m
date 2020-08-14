% Inspect the various classes identified by NMF and evaluate whether they are really coherent
% clusters

%% Set up variables

prepSR;

recs = {
    '2020-01-30/16-03-00'
    '2020-01-31/12-52-00'
    '2020-01-31/15-26-00'
    '2020-02-06/13-47-00'
    '2020-02-06/16-01-00'
    '2020-03-05/12-50-00'
    '2020-03-05/14-50-00'
    '2020-03-06/12-55-00'
    '2020-03-06/14-56-00'
    '2020-03-10/12-57-00'
    '2020-03-10/14-19-00'
    '2020-03-11/12-31-00'
    '2020-03-11/14-32-00'
    };

n_recs = length(recs);
res_mfiles = cellfun(@(r) matfile(fullfile(results_dir, r, 'mt_res.mat'), 'Writable', true), recs, 'uni', false);

n_classes = 6;
n_chans = 8;

%% Make average spectra for each class - first individual recordings

kR = randi(n_recs);
rec = recs{kR};

classes = res_mfiles{kR}.rank_classes;
freqs = res_mfiles{kR}.freq_grid;
chans = res_mfiles{kR}.name;

logz_spectra = res_mfiles{kR}.pxx_rankord;
b_nan = any(isnan(logz_spectra{1}));
logz_spectra = cellfun(@(s) s(:, ~b_nan), logz_spectra, 'uni', false);

class_spectra = zeros(length(freqs), n_classes);

figure;
tl = tiledlayout(2, 4);
title(tl, sprintf('Mean rank-order power per channel on %s', rec));

for kC = 1:n_chans
    nexttile;
    
    for kK = 1:n_classes
        class_spectra(:, kK) = mean(logz_spectra{kC}(:, classes(:, kC) == kK), 2);
    end
    
    sanePColor(1:n_classes, freqs, class_spectra, false, true);
    set(gca, 'YScale', 'log');
    xlabel('Class #');
    ylabel('Freq (Hz)');
    title(chans{kC});
end

%% Repeat above with log-z normalization

logz_spectra = res_mfiles{kR}.pxx_logz;
b_nan = any(isnan(logz_spectra{1}));
logz_spectra = cellfun(@(s) s(:, ~b_nan), logz_spectra, 'uni', false);

class_spectra = zeros(length(freqs), n_classes);

figure;
tl = tiledlayout(2, 4);
title(tl, sprintf('Mean normalized power per channel on %s', rec));

for kC = 1:n_chans
    nexttile;
    
    for kK = 1:n_classes
        class_spectra(:, kK) = mean(logz_spectra{kC}(:, classes(:, kC) == kK), 2);
    end
    
    sanePColor(1:n_classes, freqs, class_spectra, false, true);
    set(gca, 'YScale', 'log');
    xlabel('Class #');
    ylabel('Freq (Hz)');
    title(chans{kC});
end

%% Try looking at a recording with t-SNE

Y = tsne(logz_spectra{1}.');
figure;
gscatter(Y(:, 1), Y(:, 2));
