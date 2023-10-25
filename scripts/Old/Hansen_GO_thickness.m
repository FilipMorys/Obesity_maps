% This script finds the biological processes with which two gene sets are
% most involved. For each biological process (category), the mean loading
% of genes is computed (called the category score). Significance is
% assessed against a null model of category scores. This null model is
% constructed using the loadings of the same genes of interest,
% where null loadings come from PLS analysis on a
% spatial-autocorrelation preserving permutation of one matrix.

% see https://github.com/benfulcher/GeneSetEnrichmentAnalysis for more
% details on GO annotations

addpath(genpath('/dagher/dagher11/filip/Obesity_maps/software/'))
cd('/dagher/dagher11/filip/Obesity_maps/data/');

%% load
%load('gene_expression.mat') % node by gene expression matrix
%load('neurosynth.mat')      % node by term probability matrix
%load('nodes.mat')           % relevant node indices
%load('abagenes.mat')           % relevant gene indices
%load('terms.mat')           % relevant terms
load('abagenes.mat')           % all gene names
%load('spins.mat')           % spin test indices

load('../data/ProcessedData/GOAnnotationDirect-human.mat')
load('../data/ProcessedData/GOTerms_BP.mat')

%% get genes with entrezID

T = table2cell(readtable('gene_entrez_ids.csv')); % load entrezID of genes

gene_name = abagenes;     % get relevant gene names
entrezIDs = zeros(size(gene_name));

idx = [];
for k = 1:length(gene_name)                                                % for each gene
    if ismember(gene_name{k}, T(:,1))                                      % if the gene has an entrezID
        entrezIDs(k) = cell2mat(T(find(strcmp(gene_name{k}, T(:,1))),2));  % store the entrezID
        idx = [idx;k];                                                     % also store the index of the gene
    end
end

removed_genes=gene_name(entrezIDs==0);
entrezIDs = entrezIDs(entrezIDs ~= 0);                                     % remove all genes without entrezID


%% get category scores

load PLS1thickness_gene_BRs.mat %loadgene BRs from previously run PLS

PLS1geneWeightsBMIthicknessCammounLH=sortrows(PLS1geneWeightsBMIthicknessCammounLH); % sort genes alphabetically so they correspond to entrezIDs
rm_idx=cellfun(@(a) strmatch(a,PLS1geneWeightsBMIthicknessCammounLH.gene,'exact'),removed_genes,'uniform',false); %find indices of genes that were removed due to missing entrezIDs

PLS2_removedgenes=PLS1geneWeightsBMIthicknessCammounLH; %remove those genes from our PLS results
PLS2_removedgenes(cell2mat(rm_idx),:)=[];

gload=PLS2_removedgenes.BR;

ipos = find(gload > 3); % index of genes with positive loading above 3
ineg = find(gload < -3); % index of genes with negative loading below -3
gload_pos = gload(gload > 3); % loading of genes with positive loading
gload_neg = gload(gload < -3); % loading of genes with negative loading
[~,Ipos] = sort(gload_pos); % sorted
[~,Ineg] = sort(gload_neg); % sorted 

gpos_ID = entrezIDs(ipos(Ipos));   % these are the entrezIDs of the genes in the positive set
gneg_ID = entrezIDs(ineg(Ineg));  % these are the entrezIDs of the genes in the negative set

% get category scores
categoryScores = zeros(length(allGOCategories),2);
for k = 1:length(categoryScores)           % for each category
    
    GOcategory = geneEntrezAnnotations{k}; % get genes in this category
    
    category_pos = [];
    category_neg = [];
    
    for j = 1:length(GOcategory)           % for each gene in the category
        
        % check if the gene is in the positive/negative gene set
        % if so, save the loading
        if ismember(GOcategory(j),gpos_ID)
            category_pos = [category_pos; gload(find(entrezIDs == GOcategory(j)))];
        elseif ismember(GOcategory(j),gneg_ID)
            category_neg = [category_neg; gload(find(entrezIDs == GOcategory(j)))];
        end
        
    end
    
    % category score = mean loading of relevant genes in category
    categoryScores(k,1) = mean(category_pos);
    categoryScores(k,2) = mean(category_neg);   
end

%% get null category scores, run PLS again

workdir='';
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/scripts/'));
addpath(genpath('/dagher/dagher11/filip/Obesity_maps/data/'));
BMI='BMI_thickness_Cammoun_LH';
expr='Cammoun_aba_LH_thickness';
load('~/Downloads/coords125.mat');
coords=coords125(116:226,:);
bootnum=2000;
nperm=1000;
categoryScores_null = zeros(length(allGOCategories),2,nperm);

for m = 1:nperm % for each permutation
    
    % permuted PLS analysis
   
    [gload] = PLS_bootstrap_perm(expr, BMI, 5, bootnum, workdir);     % these are all BRs from a specific permutation  
    gload(cell2mat(rm_idx),:)=[];

    for k = 1:length(allGOCategories)          % for each gene in category
        
        GOcategory = geneEntrezAnnotations{k}; % get genes in the category
        
        category_pos = [];
        category_neg = [];
        
        for j = 1:length(GOcategory)           % for each gene in the category
            
            % check if the gene is in the positive/negative gene set
            % if so, save the loading
            if ismember(GOcategory(j),gpos_ID)
                category_pos = [category_pos; gload.Var2(find(entrezIDs == GOcategory(j)))];
            elseif ismember(GOcategory(j),gneg_ID)
                category_neg = [category_neg; gload.Var2(find(entrezIDs == GOcategory(j)))];
            end
        end
        
        % category score = mean loading of relevant genes in category        
        categoryScores_null(k,1,m) = mean(category_pos);
        categoryScores_null(k,2,m) = mean(category_neg);
    end
end

%% get biological processes

% get categories with a score
pos_cat = find(~isnan(categoryScores(:,1)));
neg_cat = find(~isnan(categoryScores(:,2)));

% set up p-value template
pos_pvals = zeros(length(pos_cat),1);
neg_pvals = zeros(length(neg_cat),1);

% calculate permuted p-value
for k = 1:length(pos_pvals)
    pos_pvals(k) = (1+(nnz(find(categoryScores_null(pos_cat(k),1,:) > categoryScores(pos_cat(k),1)))))/(nperm+1);
end
for k = 1:length(neg_pvals)
    neg_pvals(k) = (1+(nnz(find(categoryScores_null(neg_cat(k),2,:) < categoryScores(neg_cat(k),2)))))/(nperm+1);
end

% get biological processeses
GO_categories = table2cell(GOTable(:,2));    % these are all the category names
GO_categoriesID = table2array(GOTable(:,3)); % this is the ID of the category

positive_terms = cell(length(pos_cat),1);
negative_terms = cell(length(neg_cat),1);

% get biological process category if it exists and is available

% for the positive gene set:
for k = 1:length(pos_cat)
    positive_terms{k} = GO_categories(GO_categoriesID == allGOCategories(pos_cat(k)));
    if isempty(positive_terms{k})
        disp(allGOCategories(pos_cat(k)))
    end
end
p_idx = find(~cellfun('isempty',positive_terms));
%positive_terms = positive_terms(p_idx);           % remove empty indices
positive_terms(pos_pvals<0.05)

% for the negative gene set:
for k = 1:length(neg_cat)
    negative_terms{k} = GO_categories(GO_categoriesID == allGOCategories(neg_cat(k)));
end
n_idx = find(~cellfun('isempty',negative_terms));
%negative_terms = negative_terms(n_idx);           % remove empty indices
negative_terms(neg_pvals<0.05)

pos=table(positive_terms, pos_pvals);
pos=pos(~cellfun(@isempty, pos.positive_terms(:,1)), :);

neg=table(negative_terms, neg_pvals);
neg=neg(~cellfun(@isempty, neg.negative_terms(:,1)), :);


%% Qualitattively analyse data
figure()
pos_filtered=pos.positive_terms(pos.pos_pvals<0.001);
a=cellfun(@(x) split(x), pos_filtered, 'UniformOutput', false);
NewCellArray = vertcat( a{:} );
A = count(NewCellArray,unique(NewCellArray));
[GC, GR] = groupcounts(NewCellArray);
Word_cloud=table(GC, GR);


for i=1:height(Word_cloud)
    if length(Word_cloud.GR{i})<3
        Word_cloud.GC(i)=0;
    end
    if Word_cloud.GC(i)<10
        Word_cloud.GC(i)=0;
    end
end

wordcloud(Word_cloud, 'GR', 'GC')
exportgraphics(gcf,'Alyssas_pos_cloud.tif','Resolution',300)


figure(3)
neg_filtered=neg.negative_terms(neg.neg_pvals<0.001);
a=cellfun(@(x) split(x), neg_filtered, 'UniformOutput', false);
NewCellArray = vertcat( a{:} );
A = count(NewCellArray,unique(NewCellArray));
[GC, GR] = groupcounts(NewCellArray);
Word_cloud=table(GC, GR);

for i=1:height(Word_cloud)
    if length(Word_cloud.GR{i})<3
        Word_cloud.GC(i)=0;
    end
    if Word_cloud.GC(i)<10
        Word_cloud.GC(i)=0;
    end
end

wordcloud(Word_cloud, 'GR', 'GC')
exportgraphics(gcf,'Alyssas_neg_cloud.tif','Resolution',300)





