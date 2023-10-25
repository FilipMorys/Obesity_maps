addpath(genpath('/dagher/dagher11/filip/Software/BrainSpace-0.1.2'));
addpath(genpath('/dagher/dagher11/filip/Software/gifti-main'));
addpath(genpath('/dagher/dagher11/filip/Software/netneurotools'));
addpath(genpath('/dagher/dagher11/filip/Software/neuromaps'));


[v_dkt_lh, l_dkt_lh, ct_dkt_lh]=read_annotation('/home/bic/fmorys/nnt-data/atl-cammoun2012/fsaverage/atl-Cammoun2012_space-fsaverage_res-125_hemi-L_deterministic.annot');

%% Import data from text file
load Genes_Yengo_LH.mat

Yengo_mean=mean(genesYengoLH,2);

%% Clear temporary variables
clear opts

%% Merge labels in .annot with data file manually 

PLS2_plot=table(Yengo_mean(1:111),ct_dkt_lh.table([2:4 6:113],5));
PLS2_plot.Properties.VariableNames{1} = 'score';
PLS2_plot.Properties.VariableNames{2} = 'number';

%% Deparcellate

vertex_img=[];
unique_values=unique(l_dkt_lh);
%unique_values(6)=[]; % remove 'unknown'

for i=1:length(unique(l_dkt_lh))
    if i==18 | i==25 % CC and unknown
        continue
    end
    lab=unique_values(i);
    idx=find(l_dkt_lh==lab);
    vertex_img(idx)=PLS2_plot.score(PLS2_plot.number==lab);
end

vertex_img(vertex_img==0)=nan;

%% Plot data

plsroi=plot_hemispheres(vertex_img', {'/dagher/dagher11/filip/Obesity_maps/data/lh.white'}, 'views','lmap','labeltext',{'BMI genes'});
hold on
colormap parula %redblue

%exportgraphics(gcf,'Component2.tif','Resolution',300)