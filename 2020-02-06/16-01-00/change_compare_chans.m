% Do PCA separately for V1 and motor, then get the change speed in each
% using a small # of components and quantify order of change peaks.

%% Setup

redo_all = false; % change to true to recompute from scratch

redo_pca = redo_all || false; % change to true to redo PCA
redo_comps = redo_all || false; % change to true to re-visualize PCA and select components
redo_change = redo_all || false; % change to true to recalculate change speeds

prepSR;

date = '2020-02-06';
time = '16-01-00';
res_file = fullfile(results_dir, date, time, 'mt_res.mat');

regions = {'vis', 'mot'};

%% Do PCA, just taking 10 components b/c going to decide which ones to keep visually

pc_data = cell(2, 1);

res_mfile = matfile(res_file, 'Writable', true);

total_comps = 10;

for kR = 1:2
    res_name = sprintf('pxx_pca_1rec_%s', regions{kR});
    
    if ~redo_pca && isprop(res_mfile, res_name) % just load if possible
        pc_data(kR) = res_mfile.(res_name);
    else
        pca_opts = struct;
        pca_opts.name = sprintf('pxx_pca_1rec_%s', regions{kR});
        pca_opts.chans = kR;
        pca_opts.thresh_type = 'comps';
        pca_opts.thresh = total_comps;

        [pc_data(kR), fh] = mt_pca(res_file, pca_opts);

        savefig(fh, sprintf('pca_%s.fig', regions{kR}));
    end
end

%% Visualize to identify which components to keep (start with scatter 1-3, can do more if necessary)

comp_fname = 'change_compare_comps';

if ~redo_comps && isprop(res_mfile, comp_fname)
    comps2use = res_mfile.(comp_fname);
else
    % these times of interest have been painstakingly manually identified
    hstarts = [1521, 2197, 3722, 4530, 5180, 6081, 7619]; % in seconds
    hends = [1614, 2634, 3814, 4702, 5259, 6240, 7712];   % ditto
    
    Fs_mt = 10; % since our window step is 0.1 seconds, sample freq = 10 Hz
    htimes = arrayfun(@(s, e) Fs_mt*s+1:Fs_mt*e, hstarts, hends, 'uni', false);
    htimes = cell2mat(htimes);
    
    
    comps2use = false(2, total_comps);
    
    figure;
    
    for kR = 1:2
        comps2show = 1:3;
        
        while ~isempty(comps2show)
            
            hold off;
            dscatter(pc_data{kR}(comps2show, :));
            set(gca, 'Interactions', [zoomInteraction, rotateInteraction, rulerPanInteraction]);
            hold on;
            xlabel(sprintf('PC%d', comps2show(1)));
            ylabel(sprintf('PC%d', comps2show(2)));
            
            % show times of interest
            hdata = arrayfun(@(c) pc_data{kR}(c, htimes), comps2show, 'uni', false);
            if length(comps2show) == 3
                zlabel(sprintf('PC%d', comps2show(3)));
                scatter3(hdata{:}, 20, 'r', 'o');
            else
                scatter(hdata{:}, 20, 'r', 'o');
            end
            
            title(sprintf('Principal components in region: %s', regions{kR}));
            
            comps2add = input('Enter array of most informative components: ');
            comps2use(kR, comps2add) = true;
            
            comps2show = input('Enter 2 or 3 components to show next, or nothing to move on: ');
            while ~ismember(length(comps2show), [0, 2, 3]) || any(~ismember(comps2show, 1:total_comps))
                comps2show = input('Enter an array of 2 or 3 components, or press <enter>: ');
            end
        end        
    end
    
    % convert to 2x1 cell
    comps2use = arrayfun(@(kR) find(comps2use(kR, :)), [1; 2], 'uni', false);
    res_mfile.(comp_fname) = comps2use;
end


%% Compute change speed for each channel using pca_change

change_speed = cell(2, 1);
change_fhs = gobjects(2, 1); % figure handles

for kR = 1:2
    change_fname = sprintf('pca_change_%s', regions{kR});
    change_figname = fullfile(results_dir, date, time, [change_fname, '.fig']);
    
    if ~redo_change && isprop(res_mfile, change_fname) && exist(change_figname, 'file')
        change_speed{kR} = res_mfile.(change_fname);
        change_time = res_mfile.([change_fname, '_time']);
        change_fhs(kR) = openfig(change_figname);
        continue;
    end
    
    change_opts = struct;
    change_opts.pca_name = sprintf('pxx_pca_1rec_%s', regions{kR});
    change_opts.comps = comps2use{kR};
    change_opts.smooth_method = 'gaussian';
    change_opts.name = change_fname;
    change_opts.savefigs = false;
    
    [change_speed{kR}, change_time] = pca_change(res_file, change_opts);
    change_fhs(kR) = gcf;
    savefig(change_fhs(kR), change_figname);
end

%% Find (putative?) change peaks

ispeak = cell(2, 1);
for kR = 1:2
    data_range = max(change_speed{kR}) - min(change_speed{kR});
    ispeak{kR} = islocalmax(change_speed{kR}, 'MinProminence', 0.6 * data_range);
    
    figure(change_fhs(kR));
    hold on
    plot(change_time(ispeak{kR}), change_speed{kR}(ispeak{kR}), 'mo', 'MarkerSize', 15);
end