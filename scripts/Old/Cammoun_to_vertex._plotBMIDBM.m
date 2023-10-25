addpath(genpath('/dagher/dagher11/filip/Software/BrainSpace-0.1.2'));
addpath(genpath('/dagher/dagher11/filip/Software/gifti-main'));
addpath(genpath('/dagher/dagher11/filip/Software/netneurotools'));
addpath(genpath('/dagher/dagher11/filip/Software/neuromaps'));


[v_dkt_lh, l_dkt_lh, ct_dkt_lh]=read_annotation('/home/bic/fmorys/nnt-data/atl-cammoun2012/fsaverage/atl-Cammoun2012_space-fsaverage_res-125_hemi-L_deterministic.annot');

load BMI_DBM_LH_Cammoun.mat

BMIDBMCammoun125=BMIDBMCammoun125[1:111];

PLS2_plot=table(BMIDBMCammoun125,ct_dkt_lh.table([2:4 6:113],5));
PLS2_plot.Properties.VariableNames{1} = 'score';
PLS2_plot.Properties.VariableNames{2} = 'number';

%% Deparcellate

vertex_img=[];
unique_values=unique(l_dkt_lh);
unique_values(6)=[]; % remove 'unknown'

for i=1:length(unique(l_dkt_lh))-1
    if i==17 | i==24 % CC and unknown
        continue
    end
    lab=unique_values(i);
    idx=find(l_dkt_lh==lab);
    vertex_img(idx)=PLS2_plot.score(PLS2_plot.number==lab);
end


%% Plot data

plsroi=plot_hemispheres(vertex_img', {'/dagher/dagher11/filip/Obesity_maps/data/lh.white'}, 'views','lmap','labeltext',{'Component 1 loadings'});
hold on
%colormap redblue