% Prereq: have change_peaks from focused_change_analysis

nextpeakdelay = @(this, other) arrayfun(@(pt) min([other(other >= pt) - pt, nan]), this);
lastpeakdelay = @(this, other) arrayfun(@(pt) max([other(other <= pt) - pt, nan]), this);

figure;

regions = {'VC', 'MC'};

for kR = 1:2
    subplot(1, 2, kR); hold on;
    
    h1 = histogram(lastpeakdelay(change_peaks.real{kR}, change_peaks.real{3-kR}) / 100, 'FaceColor', 'r');
    % make sure right half has same bin width
    bwidth = h1.BinWidth;
    histogram(nextpeakdelay(change_peaks.real{kR}, change_peaks.real{3-kR}) / 100, 'FaceColor', 'b', 'BinWidth', bwidth);

    title(sprintf('Location of closest past and future %s change peaks rel to %s change peak', regions{3-kR}, regions{kR}));
    xlabel('Relative time (s)');
end