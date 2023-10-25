%% Run PLS analysis %%

%% Read in data and functions

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data/'));
BMI='BMI_DBM_LH_Cammoun';
expr='Genes_Yengo_LH';

%% Read coord_data

load('~/Downloads/coords125.mat');
coords=coords125(116:234,:);

%% Run normal PLS

bootnum=200;

PLS_bootstrap_Yengo(expr, BMI, 5, bootnum, workdir);

nperm=10000;

PLS_summary_component_stats_Yengo(expr, BMI, 5, nperm, 'spatial', workdir, coords)


