%% Run PLS analysis %%

%% Read in data and functions

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data/'));
BMI='ABCD_CT_DKT_LH';
expr='ABCDbrainspanall'; % aba_limited_expression for UKBB and HCPA, ABCD_brainspan for ABCD, HCP_brainspan for HCP

%% Read coord_data

load('DKT_coords_LH.mat');
coords=DKT_coords_LH([30,25,29,7,22,24,10,20,28,19,16],5:7);

%% Run normal PLS

bootnum=2000;

PLS_bootstrap_DKT_11regions(expr, BMI, 5, bootnum, workdir)

nperm=1000;

PLS_summary_component_stats_DKT_11regions(expr, BMI, 5, nperm, 'spatial', workdir, coords)


