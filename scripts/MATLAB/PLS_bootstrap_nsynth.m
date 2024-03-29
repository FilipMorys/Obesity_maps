function PLS_bootstrap_nsynth(GENEdata_root, MRIdata_root, ncomp, bootnum, working_dir)
% PLS regression and boostrapping
%% ------------------------- SCRIPT INFORMATION ---------------------------
%
% PLS Regression and bootstrapping
%
% This script will perform a PLS regression analysis for relative gene
% expression in 180 parcels (defined by the Glasser atlas in MNI space)
% against mean QSM signal in that group.
%
% This script is an adapted version of a script available from the
% following repository:
% https://github.com/KirstieJane/NSPN_WhitakerVertes_PNAS2016/
% Reference: Whitaker et al., 2016, PNAS
% Credit and thanks to the original author Petra Vertes
%
% Script last edited by:
% George Thomas July 2020
% - to fit with our data
% - to run with any number of components 
%
%----------------------------- INPUTS -------------------------------------
%
% GENEdata_root:
%   Name of genetic data .mat file, including Ngenes*NROI parcel expression
%   and probe information
%
% MRIdata_root:
%   Name of MRI data Nsubjects*NROI .mat file, containg participant QSM  
%   means Z-score transformed to control means for each ROI
%
% working_dir:
%   directory where MRIdata_root is and where output .csv files will go
%
% ncomp:
%   number of dimensions (components) for the PLS 
%
% bootnum:
%   number of iterations for the boostrap runs
%
%------------------------------- OUTPUTS ----------------------------------
%
% PLS<comp>_ROIscores_<MRIdata_root>.csv
%   ncomp x files containing PLS weight in NROI regions for that component
%
% PLS<comp>_geneWeights_<MRIdata_root>.csv
%   ncomp x files conatining ranked, boostrapped gene weights for each
%   component in the format 
%   geneID      geneIndex       geneZscore
%

%% ------------------------ SCIPT BEGINS HERE -----------------------------
tic
%% Import and tidy variables
disp('>>> importing + tidying variables')
disp(' ')

% get output directroy full path and add
working_dir = [working_dir];
addpath(working_dir);

% load NROI*Ngenes gene expression matrix (LH only) as well as probeInfo
%load();
GENEdata = importdata([GENEdata_root '.mat']);
%GENEdata(:,1) = [];
% replace NaN gene*ROI entries with mean expression in that ROI
%for g = 1:size(GENEdata,2)
%    GENEdata(isnan(GENEdata(:,g)),g) =...
%        mean(GENEdata(~isnan(GENEdata(:,g)),g));
%end

load('/dagher/dagher11/filip/Downloads/terms.mat');
GENEindex=(1:width(GENEdata))';
GENEids = terms.names;

% import matrix of Nsubjects*NROI MRI data
QSMdata  = importdata([MRIdata_root '.mat']);

disp('>>> grouping + averaging MRI data')
disp(' ')
% calculate average QSM score in each region and convert to column vector
%mean_MRIdata = (mean(QSMdata))';
mean_MRIdata = QSMdata;

%% run initial PLS

% DO PLS in ncomp dimensions
X=GENEdata;
Y=mean_MRIdata;
X=zscore(GENEdata,0,1);
Y=zscore(Y);
disp(['>>> running initial PLS in ' num2str(ncomp) ' dimensions'])
[~,~,XS,YS,~,PCTVAR,~,stats]=plsregress(X,Y,ncomp);
disp(YS);
disp(' ')
disp('cumulative % variance explained in Y:')
disp(' ')
disp(cumsum(100*PCTVAR(2,:)))

% calculate correlation with input data
disp('correlation of PLS components with MRI data:');
[rho, pval] = corr(mean_MRIdata, XS)

%% Plot brain values against term scores
%% Plot brain values against term scores
fig=figure('MenuBar','none','Position', [10 10 900 600]);
mdl = fitlm(XS(:,1),Y);
x=plot(mdl,'Marker','.','color','#048ba8','MarkerSize',12,'LineWidth',1,'MarkerSize',35);
%fig.WindowState = 'fullscreen'
title('')
legend('off')
xlabel('Term scores') 
ylabel('BMI*CT') 
ax = gca(fig);
ax.FontSize = 25; 
ax.Box='off';
xticks([min(xticks):0.1:max(xticks)])
ax.LineWidth=2;
set(x(2:4),'Color','#f18f01','LineWidth',4);
exportgraphics(gcf, 'Terms_nsynth_corr_new.tif','Resolution',300);



% make empty matrices to store weights and info from initial PLS run
PLSinitW = zeros(length(GENEindex),ncomp);
PLSinitIdx = zeros(length(GENEindex),ncomp);
PLSinitID = cell(length(GENEindex),ncomp);
PLSgeneIndex = zeros(length(GENEindex),ncomp);

% store info for each component and align with desired direction
for cc = 1:ncomp
    if rho(cc)<0
        stats.W(:,cc)=-1*stats.W(:,cc);
        XS(:,cc)=-1*XS(:,cc);
    end
    [PLSinitW(:,cc),PLSinitIdx(:,cc)] = sort(stats.W(:,cc),'descend');
    PLSinitID(:,cc)=GENEids(PLSinitIdx(:,cc));
    PLSgeneIndex(:,cc)=GENEindex(PLSinitIdx(:,cc));
end

%% save ROI weights to .csv
disp('>>> saving ROI weights')
for cc = 1:ncomp
    fname = [sprintf('PLS%01d_ROIscores_',cc) MRIdata_root '.csv'];
    csvwrite(fullfile(working_dir,fname),XS(:,cc));
end
disp(' ')

%% bootstrapping

% define variable for storing the (ordered) weights from all bootstrap runs
PLSbootW=zeros(length(GENEindex),bootnum,ncomp);

% start bootstrap
disp('>>> Bootstrapping - could take a while')
for bb=1:bootnum
    myresample = randsample(size(X,1),size(X,1),1);
    % define X for resampled regions
    Xr=X(myresample,:);
    % define Y for resampled regions
    Yr=Y(myresample,:);
    % perform PLS for resampled data
    [~,~,~,~,~,~,~,stats]=plsregress(Xr,Yr,ncomp);
    
    % for each component
    for cc = 1:ncomp
        % extract PLS components weights
        temp=stats.W(:,cc);
        % order the newly obtained weights the same way as initial PLS
        newW=temp(PLSinitIdx(:,cc));
        % the sign of PLS components is arbitrary - make sure this aligns between runs
        if corr(PLSinitW(:,cc),newW)<0
            newW=-1*newW;
        end
        % store (ordered) weights from this bootstrap run
        PLSbootW(:,bb,cc)=newW;
    end
    
    % report progress
    if mod(bb,100) ==0
        disp(['iteration ' num2str(bb) ' of ' num2str(bootnum)])
    end
end
disp(' ')
%% save results

% make matrix for storing stdevs and Z scores from bootstrapping
PLSstdW = zeros(length(GENEindex),ncomp);
PLSzW = zeros(length(GENEindex),ncomp);
PLSzIdx = zeros(length(GENEindex),ncomp);
PLSzID = cell(length(GENEindex),ncomp);

disp('>>> saving bootstrapped gene weights') 
% for all components
for cc = 1:ncomp
    
    % get standard deviation of weights from bootstrap runs
    PLSstdW(:,cc)=std(PLSbootW(:,:,cc),0,2);
    
    % get bootstrap weights
    temp=PLSinitW(:,cc)./PLSstdW(:,cc);
    
    % order bootstrap weights (Z) and names of genes
    [PLSzW(:,cc), PLSzIdx(:,cc)]=sort(temp,'descend');
    PLSzID(:,cc)=PLSinitID(PLSzIdx(:,cc));
    PLSgeneIndex(:,cc)=PLSgeneIndex(PLSzIdx(:,cc));
    
    % save results to .csv
    fname = [sprintf('PLS%01d_geneWeights_',cc) MRIdata_root '.csv'];
    fid = fopen(fullfile(working_dir,fname),'w');
    for ii=1:length(GENEids)
        fprintf(fid,'%s, %d, %f\n',...
            PLSzID{ii,cc}, PLSgeneIndex(ii,cc), PLSzW(ii,cc));
    end
    fclose(fid);
end
disp(' ')
disp('Done!')
%% report speed
toc