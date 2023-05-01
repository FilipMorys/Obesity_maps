addpath(genpath('/dagher/dagher11/filip/Software/BrainSpace-0.1.2'));
addpath(genpath('/dagher/dagher11/filip/Software/gifti-main'));
addpath(genpath('/dagher/dagher11/filip/Software/netneurotools'));
addpath(genpath('/dagher/dagher11/filip/Software/neuromaps'));


[v_dkt_lh, l_dkt_lh, ct_dkt_lh]=read_annotation('/dagher/dagher11/filip/Obesity_maps/data_old/DKT.annot');
[v_dkt_rh, l_dkt_rh, ct_dkt_rh]=read_annotation('/dagher/dagher11/filip/Obesity_maps/data_old/rh.DKT.annot');


load("HCPA_CT_DKT_RH.mat");
load("HCPA_CT_DKT_LH.mat");


%% Clear temporary variables
clear opts

%% Merge labels in .annot with data file manually 

PLS1_plot=table(HCPA_CT_DKT_LH',ct_dkt_lh.table([3:4,6:32, 35:36],5));
PLS1_plot.Properties.VariableNames{1} = 'score';
PLS1_plot.Properties.VariableNames{2} = 'number';


PLS1_plotrh=table(HCPA_CT_DKT_RH',ct_dkt_lh.table([3:4,6:32, 35:36],5));
PLS1_plotrh.Properties.VariableNames{1} = 'score';
PLS1_plotrh.Properties.VariableNames{2} = 'number';

%% Deparcellate

vertex_img=[];
unique_values=unique(l_dkt_lh);
unique_values(6)=[]; % remove 'unknown'

for i=1:length(unique(l_dkt_lh))-1
    lab=unique_values(i);
    idx=find(l_dkt_lh==lab);
    vertex_img(idx)=PLS1_plot.score(PLS1_plot.number==lab);
end

vertex_imgrh=[];
unique_values=unique(l_dkt_rh);
unique_values(6)=[]; % remove 'unknown'

for i=1:length(unique(l_dkt_rh))-1
    lab=unique_values(i);
    idx=find(l_dkt_rh==lab);
    vertex_imgrh(idx)=PLS1_plotrh.score(PLS1_plotrh.number==lab);
end

%% Plot data

plsroi=plot_hemispheres([vertex_img'; vertex_imgrh'], {'/dagher/dagher11/filip/Obesity_maps/data_old/lh.white','/dagher/dagher11/filip/Obesity_maps/data_old/rh.white'}, 'views','lmap','labeltext',{'BMI*CT'});
%plsroi.colorlimits([-3,3]);
hold on
colormap redblue
%for ABCD
cmap=colormap(redblue);
colormap(cmap(1:128,:));

exportgraphics(gcf,'/dagher/dagher11/filip/Obesity_maps/HCPA/BMI.tif','Resolution',300)