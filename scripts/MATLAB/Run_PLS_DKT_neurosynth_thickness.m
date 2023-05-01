%% Run PLS analysis %%

%% Read in data and functions

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data_old/'));

BMI='ABCD_CT_DKT_LHRH';
nsynth='Nsynth_DKT_thickness_LHRH';

%% Read coord_data

load('DKT_coords_LHRH.mat');
coords=DKT_coords(:,5:7);

%% Run normal PLS

bootnum=20000;

PLS_bootstrap_nsynth(nsynth, BMI, 5, bootnum, workdir);

nperm=10000;

PLS_summary_component_stats_nsynth_thickness(nsynth, BMI, 5, nperm, 'spatial', workdir, coords)


