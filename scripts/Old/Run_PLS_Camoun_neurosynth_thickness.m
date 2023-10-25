%% Run PLS analysis %%

%% Read in data and functions

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data/'));
BMI='BMI_Cammoun_thickness_LHRH';
nsynth='Nsynth_Cammoun_thickness_LHRH';

%% Read coord_data

load('~/Downloads/coords125.mat');
coords=coords125([1:108,116:226],:);

%% Run normal PLS

bootnum=20000;

PLS_bootstrap_nsynth(nsynth, BMI, 5, bootnum, workdir);

nperm=10000;

PLS_summary_component_stats_nsynth_thickness(nsynth, BMI, 5, nperm, 'spatial', workdir, coords)


