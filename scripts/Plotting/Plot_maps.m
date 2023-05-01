addpath(genpath('/dagher/dagher11/filip/Software/BrainSpace-0.1.2'));
addpath(genpath('/dagher/dagher11/filip/Software/gifti-main'));
addpath(genpath('/dagher/dagher11/filip/Software/netneurotools'));
addpath(genpath('/dagher/dagher11/filip/Software/neuromaps'));
addpath(genpath('/dagher/dagher11/filip/Software/nifti_toolbox'));

close all
BMI_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/BMI_rh.gii');
BMI_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/BMI_lh.gii');

dasb_2017_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/beliveau2017dasb_L.gii');
dasb_2017_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/beliveau2017dasb_R.gii');

dasb_2012_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/savli2012dasb_L.gii');
dasb_2012_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/savli2012dasb_R.gii');

dasb_2012p_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/savli2012p943_L.gii');
dasb_2012p_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/savli2012p943_R.gii');

fpcit_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/dukart2018fpcit_L.gii');
fpcit_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/dukart2018fpcit_R.gii');

ec1r_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/laurikainen2018fmpepd2_L.gii');
ec1r_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/laurikainen2018fmpepd2_R.gii');

BMI_area_rh=gifti('/dagher/dagher11/filip/Obesity_maps/data/BMI_area_rh.gii');
BMI_area_lh=gifti('/dagher/dagher11/filip/Obesity_maps/data/BMI_area_lh.gii');

BMI_both=zscore([BMI_lh.cdata; BMI_rh.cdata]);
dasb_both=zscore([dasb_2017_lh.cdata; dasb_2017_rh.cdata]);
dasb2012_both=zscore([dasb_2012_lh.cdata; dasb_2012_rh.cdata]);
dasb2012p_both=zscore([dasb_2012p_lh.cdata; dasb_2012p_rh.cdata]);
fpcit_both=zscore([fpcit_lh.cdata; fpcit_rh.cdata]);
ec1r_both=zscore([ec1r_lh.cdata; ec1r_rh.cdata]);
BMIarea_both=zscore([BMI_area_lh.cdata; BMI_area_rh.cdata]);


plot=[BMI_both, dasb_both, dasb2012_both, fpcit_both, ec1r_both];
plot_hemispheres(plot, {'/dagher/dagher11/filip/Obesity_maps/data/lh.white', '/dagher/dagher11/filip/Obesity_maps/data/rh.white'}, 'views','lmap','labeltext',{'BMI - CT','5-HTT','p943','DAT','CB1R'});

hold on
%colormap redblue %redblue

exportgraphics(gcf,'Brain_maps_thickness.tif','Resolution',300)

%% Plot gene map for LH from Yengo 2018

%genes=nifti('/home/bic/fmorys/Desktop/image.nii');

