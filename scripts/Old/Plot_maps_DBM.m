addpath(genpath('/dagher/dagher11/filip/Software/BrainSpace-0.1.2'));
addpath(genpath('/dagher/dagher11/filip/Software/gifti-main'));
addpath(genpath('/dagher/dagher11/filip/Software/netneurotools'));
addpath(genpath('/dagher/dagher11/filip/Software/neuromaps'));
addpath(genpath('/dagher/dagher11/filip/Software/nifti_toolbox'));

close all

load BMI_DBM.mat

dasb_2017_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/beliveau2017dasb_L.gii');
dasb_2017_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/beliveau2017dasb_R.gii');

ec1r_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/laurikainen2018fmpepd2_L.gii');
ec1r_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/laurikainen2018fmpepd2_R.gii');

flubatine_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/hillmer2016flubatine_L.gii');
flubatine_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/hillmer2016flubatine_R.gii');

p943_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/savli2012p943_L.gii');
p943_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/savli2012p943_R.gii');

cmrglu_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/raichlecmruglu_L.gii');
cmrglu_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/raichlecmruglu_R.gii');

BMI_both=zscore([BMIDBMLH.e04; BMIDBMRH.e04]);
dasb_both=zscore([dasb_2017_lh.cdata; dasb_2017_rh.cdata]);
ec1r_both=zscore([ec1r_lh.cdata; ec1r_lh.cdata]);
flubatine_both=zscore([flubatine_lh.cdata; flubatine_rh.cdata]);
p943_both=zscore([p943_lh.cdata; p943_rh.cdata]);
cmrglu_both=zscore([cmrglu_lh.cdata; cmrglu_rh.cdata]);

plot=[BMI_both];
plot_hemispheres(plot, {'/dagher/dagher11/filip/Obesity_maps/data/lh.white', '/dagher/dagher11/filip/Obesity_maps/data/rh.white'}, 'views','lmap','labeltext',{'BMI - DBM'});

hold on
colormap parula %redblue

exportgraphics(gcf,'DBM_map.tif','Resolution',300)

figure()
plot=[dasb_both, ec1r_both, flubatine_both, p943_both, cmrglu_both];
plot_hemispheres(plot, {'/dagher/dagher11/filip/Obesity_maps/data/lh.white', '/dagher/dagher11/filip/Obesity_maps/data/rh.white'}, 'views','lmap','labeltext',{'5-HTT','CB1R','nAChR','5-HT1B','CMR Glc'});

hold on
colormap parula %redblue

exportgraphics(gcf,'Brain_maps.tif','Resolution',300)
%% Plot gene map for LH from Yengo 2018

%genes=nifti('/home/bic/fmorys/Desktop/image.nii');