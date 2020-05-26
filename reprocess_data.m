%% 2020-01-31/15-26-00

prepSR;

recdate = '2020-01-31';
time = '15-26-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;

options.chans = [11, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'};

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

options.save = false;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-02-06/13-47-00

prepSR;

recdate = '2020-02-06';
time = '13-47-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.artifacts = [
    1387, 1391
    5050, len_secs % (to end of recording)
    ];

options.chans = [11, 22, 43, 55];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-02-06/16-01-00

prepSR;

recdate = '2020-02-06';
time = '16-01-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.savedir = savedir;
options.artifacts = [
    832, 834
    ];

options.chans = [9, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'};

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-05/12-50-00

prepSR;

recdate = '2020-03-05';
time = '12-50-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [];

options.chans = [11, 22, 43, 55];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-05/14-50-00

prepSR;

recdate = '2020-03-05';
time = '14-50-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [];

options.chans = [11, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-06/12-55-00

prepSR;

recdate = '2020-03-06';
time = '12-55-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    4849, 4859      % "black swan" 1
    4889, 4899      % "black swan" 2
    ];

options.chans = [11, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-06/14-56-00

prepSR;

recdate = '2020-03-06';
time = '14-56-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    5454, 5455
    6324, 6326
    ];

options.chans = [11, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'};                   
                  
options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-10/12-57-00

prepSR;

recdate = '2020-03-10';
time = '12-57-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [];

options.chans = [10, 22, 43, 55];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-10/14-19-00

prepSR;

recdate = '2020-03-10';
time = '14-19-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    8673.75, 8675
];

options.chans = [11, 23, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-11/12-31-00

prepSR;

recdate = '2020-03-11';
time = '12-31-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    288, 290
];

options.chans = [11, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'}; 

options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);

%% 2020-03-11/14-32-00

prepSR;

recdate = '2020-03-11';
time = '14-32-00';

savedir = fullfile(results_dir, recdate, time);

data_s = load(fullfile(processed_lfp_dir, sprintf('meanSub_%s_%s.mat', recdate, time)));

len_secs = size(data_s.meanSubFullTrace, 2) / data_s.finalSampR;

options = struct;
options.artifacts = [
    3039, 3041
    4120, 4169
    4194, 4211
    4222, 4225
    4736, 4738
    6109, 6111
];

options.chans = [11, 22, 43, 54];
options.chan_names = {'V1 mid-super', 'V1 mid-deep', 'MC mid-super', 'MC mid-deep'};                   
                  
options.save = false;

% smaller window
options.window = 6;
options.padbase = 60; % use padding as if it were 60 seconds
options.winstep = 0.1;

mt_res = multitaper_analysis(data_s, options);

chan_fixup(mt_res, savedir);
