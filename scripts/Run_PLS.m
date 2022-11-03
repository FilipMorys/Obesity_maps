%% Run PLS analysis %%

%% Read in data and functions

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
BMI='BMI_DBM_LH';
expr='DKT_aba_LH';

%% Read coord_data

load('DKT_coords_LH.mat');
coords=DKT_coords_LH(:,5:7);

%% Run normal PLS

bootnum=20000;

PLS_bootstrap(expr, BMI, 5, bootnum, workdir)

nperm=10000;

PLS_summary_component_stats('abagenesCammoun125', BMI, 5, nperm, 'spatial', workdir, coords)


