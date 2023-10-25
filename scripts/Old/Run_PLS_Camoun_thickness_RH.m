%% Run PLS analysis %%

%% Read in data and functions

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data/'));
BMI='BMI_thickness_Cammoun_RH';
expr='Cammoun_aba_RH_thickness';

%% Read coord_data

load('~/Downloads/coords125.mat');
coords=coords125(1:108,:);

%% Run normal PLS

bootnum=20000;

PLS_bootstrap(expr, BMI, 5, bootnum, workdir)

nperm=10000;

PLS_summary_component_stats_thick_RH(expr, BMI, 5, nperm, 'spatial', workdir, coords)


