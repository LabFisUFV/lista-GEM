%% PROTOCOL FOR THE RECONSTRUCTION OF A Lipomyces starkeyi GENOME SCALE MODEL USING THE RAVEN TOOLBOX
%  AUTHORS: Eduardo L. M. Almeida; Maurício A. M. Ferreira
%
% Complementary to the manuscript, where more detailed descriptions are provided. 
%
% This protocol was prepared based on the reconstruction protocols and scripts developed for
% Rhodotorula toruloides and Yarrowia lipolytica  both available on GitHub:
% https://github.com/SysBioChalmers/rhto-GEM and https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM.
%
% Define paths for the folders containing the files needed or generated during
% the model reconstruction. Run this script every time you start a new MATLAB
% session to work with this script.
clear; clc;
if ~exist([pwd() '/reconstructionProtocol.m']); error(['Make sure that '...
        'your Current Folder is the one containing the reconstructionProtocol file.']); end
cd ../;  root = [pwd() '/'];
data = [root 'data/'];
code = [root 'code/'];
cd(code)

%% %% 1
% This script generates a first draft of a genome-scale model (here: L.
% starkeyi) using the getModelFromHomology function in RAVEN. It uses
% R. toruloides as a template model.

% BLAST against Rhto genome. L. starkeyi protein fasta obtained from NCBI (ID: 10576)
% blastRhto  = getBlast('lista',[data '\genomes\lista.fasta'],'rhto',[data '\genomes\rhto.fasta']);
% blastYli   = getBlast('lista',[data '\genomes\lista.fasta'],'yli',[data 'genomes\iyali.fasta']);
% 
% mkdir([root '/scrap'])
% 
% save([root '/scrap/blastStruct.mat'],'blast*');
load([root '/scrap/blastStruct.mat']);

% Load R. toruloides model, downloaded from 
% https://github.com/SysBioChalmers/rhto-GEM (version 1.3.0).
modelRhto        = importModel([data '/templateModels/rhto.xml'],true);
modelRhto.id     = 'rhto';
% Remove compartment information from metabolite ids (redundant).
modelRhto.mets   = regexprep(modelRhto.mets,'\[[a-z]+\]$','');

% Change reversibility fields to match boundaries. Prevents problems with
% MENECO.
% Rhto
for i=1:length(modelRhto.rxns)
	if modelRhto.lb(i)==0 && not(modelRhto.ub(i)==0)
        modelRhto.rev(i)=0;
    elseif modelRhto.lb(i)<=0 && modelRhto.ub(i)>=0
        modelRhto.rev(i)=1;
	end
end

% Confirm that the model is functional, set objective to growth.
modelRhto = setParam(modelRhto,'obj','r_4041',1);
solveLP(modelRhto)

% Yarrowia 
modelYli        = importModel([data '/templateModels/iYali.xml'],true);
modelYli.id     = 'yli';
[modelYli.grRules, modelYli.rxnGeneMat]=standardizeGrRules(modelYli);
for i=1:length(modelYli.rxns)
	if modelYli.lb(i)==0 && not(modelYli.ub(i)==0)
        modelYli.rev(i)=0;
    elseif modelYli.lb(i)<=0 && modelYli.ub(i)>=0
        modelYli.rev(i)=1;
	end
end
solveLP(modelYli)

modelYli.rxns = regexprep(modelYli.rxns,'y00','r_');
modelYli = removeReactions(modelYli,contains(modelYli.rxns,'y10'),true,true,true);

save([root '/scrap/modelTemplate.mat'], 'model*');

%% Generate draft model, based on homology.
model   = getModelFromHomology(modelRhto,blastRhto,'lista',{},1,false,10^-20,150,35);

%% Add some reactions as based on homology with Yarrowia lipolytica
modelYli    = getModelFromHomology(modelYli,blastYli,'lista',{},1,false,10^-20,150,35);

% Discard reactions that were already in draft lista-GEM
modelYli    = removeReactions(modelYli,contains(modelYli.rxns,model.rxns),true,true,true);

% Focus on reactions derived from rhto-GEM
tmp         = removeReactions(modelYli,cellfun(@isempty,regexp(modelYli.rxns,'r_\d{4}$')),true,true,true);

% How the Yarrowia model was constructed, there is a set of new metabolites
% that were introduced by simplifying lipid metabolism. Discard these
% metabolites and associated reactions.
tmp         = removeMets(tmp,contains(tmp.mets,'m'),false,true,true,true);
tmp         = removeReactions(tmp,~ismember(tmp.rxns,modelRhto.rxns));
model       = addRxnsGenesMets(model,modelRhto,tmp.rxns,tmp.grRules,'Identified from homology to Yarrowia lipolytica',2);

% Add Yarrowia specific reactions
tmp                 = removeReactions(modelYli,~contains(modelYli.rxns,'y'),true,true,true);
% Replace old identifiers with new format
oldIdx              = find(contains(tmp.mets,'m'));
tmp.mets(oldIdx)    = generateNewIds(model,'mets','m_',numel(oldIdx));
tmp.metNames        = regexprep(tmp.metNames,'alpha-D-ribose 1-phosphate','alpha-D-ribose 1-phosphate(2-)');

modelComb           = mergeModels({model,tmp});
model               = contractModel(modelComb);

%% Add L. starkeyi GEM meta data
model.annotation.defaultLB    = -1000;
model.annotation.defaultUB    = +1000;
model.annotation.taxonomy     = 'taxonomy/675824';
model.annotation.givenName    = 'Eduardo Almeida & Mauricio Ferreira';
model.annotation.familyName   = 'LABFIS/UFV';
model.annotation.email        = 'eduardo.menezes@ufv.br & mauricio.moura@ufv.br';
model.annotation.organization = 'Universidade Federal de Vicosa';
model.annotation.note         = 'Lipomyces starkeyi NRRL Y-11557';
model.id                      = 'lista';
model.description             = 'Genome-scale metabolic model of Lipomyces starkeyi';

save([root '/scrap/model_r1.mat'],'model');

% To inspect the first draft model:
exportToExcelFormat(model,[root '/scrap/model_r1.xlsx']);

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model); cd('reconstruction')

%% %% 2
%% Copy pseudoreactions
load([root '/scrap/model_r1.mat']);
load([root '/scrap/modelTemplate.mat']);

rxns    = modelRhto.rxns(contains(modelRhto.rxnNames,'pseudoreaction'));
model   = addRxnsGenesMets(model,modelRhto,rxns,false,'Modeling reaction',1);
model   = addRxnsGenesMets(model,modelRhto,'r_4046',false,'Modeling reaction',1);
rxns    = modelRhto.rxns(contains(modelRhto.rxnNames,'SLIME'));
model   = addRxnsGenesMets(model,modelRhto,rxns,false,'SLIME reaction',1);

%% Add all exchange rxns
% These were not gene annotated, and therefore not added in draft.
% Might not require all exchange rxns, but easier to remove unconnected ones later.
rxns    = getExchangeRxns(modelRhto);
model   = addRxnsGenesMets(model,modelRhto,rxns,false,'Modelling reaction',1);

%% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx   = find(cellfun(@isempty,modelRhto.grRules)); % Which rxns have no genes
rxnIdx      = find(getTransportRxns(modelRhto));
rxnIdx      = intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns        = modelRhto.rxns(rxnIdx); % Obtain reaction IDs
model       = addRxnsGenesMets(model,modelRhto,rxns,false,'Modeling reaction required for intercellular transport, gene unknown',1);

save([root '/scrap/model_r2.mat'],'model');

% To inspect the second draft model:
exportToExcelFormat(model,[root '/scrap/model_r2.xlsx']);

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model);cd('reconstruction')

