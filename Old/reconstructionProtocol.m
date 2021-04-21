%% PROTOCOL FOR THE RECONSTRUCTION OF A Lipomyces starkeyi GENOME SCALE MODEL USING THE RAVEN TOOLBOX
%  AUTHORS: Eduardo L. M. Almdeida; Maurício A. M. Ferreira
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


%% IMPORT TEMPLATE MODEL
% Load Y. lipolytica template GEM using importModel()
% Source: https://github.com/SysBioChalmers/Yarrowia_lipolytica_W29-GEM
modeliYali = importModel([data 'templateModels/iYali.xml'], true);
modeliYali.id = 'iYali';

% If one wants to evaluate what the model contains, it can be easier to
% browse through it as an Excel sheet. During the reconstruction we
% occasionally want to write files that we don't necessarily want to keep
% and track, but just make temporarily. We can write these to a 'scrap'
% folder that will not be tracked on GitHub.
mkdir([root 'scrap'])
exportToExcelFormat(modeliYali, [root 'scrap/iYali.xlsx']);

%{
Note that RAVEN is used to write iYali in its repository, so there is
no redundant indication of compartments in the metabolite IDs that we
otherwise would have wanted to clear out.
%}

% It can be useful to store intermediate states of the MATLAB environment.
% We will store these in the 'scrap' folder, loading them the next time
% will then allow us to continue from this point.
save([root 'scrap/importModel.mat'])
% load([root 'scrap/importModel.mat'])

%% GENERATE MODELS FROM HOMOLOGY 
% MATCH PROTEIN FASTA IDs
% DETERMINE HOMOLOGY by BLAST 

% BLAST the whole-genome protein FASTA of L. starkeyi against the
% Y. lipolytica protein FASTA. This can take some minutes.
blastiYali = getBlast('lista',[data '/genomes/lista.fasta'], ...
            'iYali',[data '/genomes/iYali.faa']);
 
% Save intermediate files in 'scrap' folder that is not tracked by Git
save([root '/scrap/blastStruct.mat'],'blast*');

% The blast results are then used to generate a new draft model.
model=getModelFromHomology(modeliYali,blastiYali,'lista',{},1,false,10^-20,150,35);

save([root '/scrap/model_r1.mat'],'model');
%load([root 'scrap/model_r1.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% To inspect the first draft model:
exportToExcelFormat(model,[root '/scrap/r1_listaGEM.xlsx']);

% Add exchange reactions for media components
mediumComps = {'r_1654', 'r_1672', 'r_1808', 'r_1832', 'r_1861', ...
               'r_1992', 'r_2005', 'r_2060', 'r_2100', 'r_2111'};
model = addRxnsGenesMets(model, modelRhto, mediumComps);

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Extport to inspect:
exportToExcelFormat(model, [root 'scrap/r2_paplaGEM.xlsx']);

% Save workspace
save([root 'scrap/homology.mat'])
% load([root 'scrap/homology.mat'])
clear blast mediumComps

%% DEFINE BIOMASS COMPOSITION
% P. laurentii has poly-unsaturated fatty acids, similar to R. toruloides.
% Use the biomass pseudoreactions from rhto-GEM as template to modify.

% Find all reactions with 'pseudreaction' in reactio name in rhto-GEM, and
% add these to the draft model.
biomassRxns = modelRhto.rxns(endsWith(modelRhto.rxnNames, 'pseudoreaction'));
model = addRxnsGenesMets(model, modelRhto, biomassRxns);

% Add all exchange rxns
% These were not gene annotated, and therefore not added in draft.
% Might not require all exchange rxns, but easier to remove unconnected ones later.
rxns    = getExchangeRxns(modelRhto);
model   = addRxnsGenesMets(model,modelRhto,rxns,false,'Modelling reaction',1);

% Same as exchange reactions, add all non-gene annotated transport reactions
noGeneIdx   = find(cellfun(@isempty,modelRhto.grRules)); % Which rxns have no genes
rxnIdx      = find(getTransportRxns(modelRhto));
rxnIdx      = intersect(rxnIdx,noGeneIdx); % Keep the ones without gene notation
rxns        = modelRhto.rxns(rxnIdx); % Obtain reaction IDs
model       = addRxnsGenesMets(model,modelRhto,rxns,false,'Modeling reaction required for intercellular transport, gene unknown',1);

% For the lipid curation and gapfilling 
model = addRxnsGenesMets(model, modelRhto,{'r_4062', 'r_4064', 'r_4046'});
model = setParam(model, 'ub', {'r_4062', 'r_4064', 'r_4046'}, 1000);
model = setParam(model, 'lb', {'r_4062', 'r_4064', 'r_4046'}, 0);

% Load biomass information
fid           = fopen([data 'biomass/biomassCuration.csv']);
loadedData    = textscan(fid, '%q %q %q %f','delimiter', ',', 'HeaderLines', 1);
fclose(fid);

BM.name       = loadedData{1};    BM.mets     = loadedData{2};
BM.pseudorxn  = loadedData{3};    BM.coeff    = loadedData{4};

% Nucleotides (DNA)
% Find out which rows contain the relevant information
indexes = find(contains(BM.pseudorxn, 'DNA'));
% Define new stoichiometries
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
% Change reaction
model = changeRxns(model, 'r_4050', equations, 1);

% Ribonucleotides (RNA)
indexes = find(contains(BM.pseudorxn, 'RNA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4049', equations, 1);

% Amino acids (protein)
indexes = find(contains(BM.pseudorxn, 'AA'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4047', equations, 1);

% Carbohydrates
indexes = find(contains(BM.pseudorxn, 'carbohydrate'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4048', equations, 1);

% Lipid backbones
indexes = find(contains(BM.pseudorxn, 'backbone'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4063', equations, 1);

% Lipid chains
indexes = find(contains(BM.pseudorxn, 'chain'));
equations.mets          = BM.mets(indexes);
equations.stoichCoeffs  = BM.coeff(indexes);
model = changeRxns(model, 'r_4065', equations, 1);

save([root 'scrap/biomass.mat'])
% load([root 'scrap/biomass.mat'])
clear indexes equations loadedData fid BM biomassRxns ans

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Extport to inspect:
exportToExcelFormat(model, [root 'scrap/r3_paplaGEM.xlsx']);


%% CURATION OF LIPID REACTIONS
% P. laurentii has unique fatty acid and lipid class compositions. SLIMEr
% explicitly models each lipid moiety, with unique chain distribution, but
% to reduce complexity we will only include a subset of possible chain
% distributions. To do this, files with templates reactions will be
% modified to match the desired chain distributions. First read the file.
fid         = fopen([data '/reconstruction/lipidTemplates.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', ...
    '\t', 'HeaderLines', 1);
fclose(fid);

% Reorganize the content so that it can be used by the addLipidReactions
% function.
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.grRules     = loadedData{4};
template.chains = {};
for k = 1:length(loadedData)-4; template.chains(:,k) = loadedData{k+4}; end


% Remove reactions that match a lipid template reaction (ignoring acyl-chains)
toRemove    = regexprep(template.rxns, 'CHAIN.*', '');
toRemove    = find(startsWith(model.rxnNames, toRemove));
model       = removeReactions(model, toRemove);

% Now use the templates to add the relevant reactions to the model. If a
% reaction already existed in the R. toruloides template model, then it
% will use the same reaction identifier.
cd([code 'lipidMetabolism'])
model = addLipidReactions(template, model, modelRhto);

fid         = fopen([data '/reconstruction/lipidTransport.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', ...
    '\t', 'HeaderLines', 1);
fclose(fid);
clear template
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+3}; end

model = addLipidReactions(template, model, modelRhto);


% Apply SLIME reactions
% First remove any SLIME reactions that might exist in the draft model.
model = removeReactions(model,contains(model.rxnNames,'SLIME rxn'));

% Load SLIME template reactions and parse it through the addSLIMEreactions
% function to amend the model.
fid             = fopen([data '/reconstruction/SLIMERtemplates.txt']);
firstLine       = fgets(fid);
numCols         = numel(strfind(firstLine,char(9))); % number of \t
loadedData      = textscan(fid,['%q %q %f' repmat(' %q',[1,numCols-2])],'delimiter','\t');
fclose(fid);
template.metName    = loadedData{1};	template.bbID   = loadedData{2};
template.bbMW       = loadedData{3};    template.comps  = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');

model = addSLIMEreactions(template, model, modelRhto);
cd(code)

save([root 'scrap/lipids.mat'])
% load([root 'scrap/lipids.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Extport to inspect:
exportToExcelFormat(model, [root 'scrap/r4_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r4_paplaGEM.xml'])

clear fid loadedData template k toRemove ans

%% 3.6 PERFORM GAP-FILLING
% Use biomass production as obj func for gapfilling
model = setParam(model, 'obj', 'r_4041', 1);

% Set biomass production to arbitrary low flux, to force gap-filling to
% produce biomass.
model = setParam(model, 'lb', 'r_4041', 0.01);

% Set glycose uptake at a higher value, to make sure that fluxes are high
% enough during gap-filling that they won't be ignored due to the tolerance
% of the MILP solver
model = setParam(model, 'lb', 'r_1714', -1);

% From the Rhto model, remove all exchange reactions (the
% necessary ones we already added, don't want to add new ones)
modelRhto2 = removeReactions(modelRhto, getExchangeRxns(modelRhto, 'both'), true, true, true);

% Run fillGaps function
[~, ~, addedRxns, model] = fillGaps(model, modelRhto2, false, true);

% Verify that model can now grow
sol = solveLP(model, 1)
printFluxes(model, sol.x)

cd([code 'lipidMetabolism'])
model = scaleLipids(model, 'tails');
cd(code)

% Verify the change of fluxes in the model after lipid scaling
sol = solveLP(model, 1)
printFluxes(model, sol.x)

% Set maximum uptake of carbon source back to 1 mmol/gDCW*hr
model = setParam(model, 'lb', 'r_1714', -1);
model = setParam(model, 'lb', 'r_4041', 0);
model = deleteUnusedGenes(model);

save([root '/scrap/gapfilling.mat'],'model');
% load([root 'scrap/gapfilling.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/r5_paplaGEM.xlsx']);

clear addedRxns rxns sol biomassRxns


%% 3.7 MANUAL CURATION
% Manual curation identified some more reactions, e.g. xylulokinase and
% complex IV were missing.

fid         = fopen([data '/reconstruction/manualCuration.txt']);
loadedData  = textscan(fid,'%q %q','delimiter','\t'); fclose(fid);
rxns        = loadedData{1};
grRules     = regexprep(loadedData{2},'***','');

model = addRxnsGenesMets(model,modelRhto,rxns,grRules,'Identified from homology, manual curation',2);


% Modify gene associations of gap-filled reactions. Find reactions
% annotationed with R. toruloides genes.
rhtoRxns = find(contains(model.grRules, 'RHTO'));
model.rxnNames(rhtoRxns);
model.grRules(rhtoRxns);

fid         = fopen([data '/reconstruction/updateGrRules.txt']);
loadedData  = textscan(fid,'%q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

rxns = loadedData{1};
grRules = loadedData{2};
model = changeGrRules(model,rxns,grRules);

    
% Remove rhto remnants in subSystems
% As a remnant of the homology based reconstruction, some of the more
% complex grRules have redundancies in subunit configurations. Also remove
% unused metabolites and remove 'rhto' prefix from subsystems.
run([code 'curation/cleanupModel']);

% Set exchange reactions to alternative carbon sources to reversible
model = setParam(model,'rev',{'r_1634','r_1706','r_1710','r_1711','r_1714','r_1715','r_1718','r_1808','r_2058'},1);

% Set ICDH, THFS and all fatty-acid-CoA ligases as irreversible
model = setParam(model,'lb',{'r_0659','r_0446'},0);
model = setParam(model,'lb',contains(model.rxnNames,'fatty-acid--CoA ligase'),0);
model = setParam(model,'ub','r_4046',1000);

% Save workspace
save([root 'scrap/cleanup.mat'])
% load([root 'scrap/cleanup.mat'])

%Extport to inspect:
exportToExcelFormat(model, [root '/scrap/r6_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/r6_paplaGEM.xml'])


disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
 


%% 3.8 SIMULATIONS
% Parameters for simulating batch growth in minimal medium
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));

rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});

targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));

% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);
%
% set glucose as carbon source and unlimited O2 for aerobic growth
model = setParam(model, 'lb', {'r_1714'}, 0);   % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);

model = setParam(model, 'lb', {'r_1634'}, 0);   % acetato
model = setParam(model, 'ub', {'r_1634'}, 0);
model = setParam(model, 'lb', {'r_1710'}, 0);   % galactose
model = setParam(model, 'ub', {'r_1710'}, 0);
model = setParam(model, 'lb', {'r_2058'} ,0);   % sacarose
model = setParam(model, 'ub', {'r_2058'}, 0);
model = setParam(model, 'lb', {'r_1718'}, 0);   % xilose
model = setParam(model, 'ub', {'r_1718'}, 0);
model = setParam(model, 'lb', {'r_1711'}, 0);   % ac galacturônico
model = setParam(model, 'ub', {'r_1711'}, 0);
model = setParam(model, 'lb', {'r_1808'}, 0);   % glicerol
model = setParam(model, 'ub', {'r_1808'}, 0);
model = setParam(model, 'lb', {'r_1706'}, 0);   % arabinose
model = setParam(model, 'ub', {'r_1706'}, 0);
model = setParam(model, 'lb', {'r_1715'}, -2.01);   % manose
model = setParam(model, 'ub', {'r_1715'}, 0);

model = setParam(model, 'lb', {'r_1992'}, -1000);   % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model = setParam(model, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model = setParam(model, 'ub', {'r_2111'}, 1000);
model = setParam(model, 'obj',{'r_2111'}, 1);   
%
sol = solveLP(model,1)
%
printFluxes(model, sol.x, true);
%
%printFluxes(model, sol.x, false);

%%
% nutrient uptake reactions to simulate complex medium conditions
modelYPD = model;

aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
             
modelYPD = setParam(modelYPD, 'lb', aminoacidRxns, -0.01);

solYPD = solveLP(modelYPD)

printFluxes(modelYPD, solYPD.x, true);

%printFluxes(modelYPD, solYPD.x, false);

%% Chemostat simulation
% adjust parameters for chemostat growth
model = setParam(model, 'eq', {'r_2111'}, 0.05);  % fix specific growth rate at the dilution rate value

uptake = find(strcmp(model.rxnNames,'D-glucose exchange')); % remove constraints on substrate uptake
model = setParam(model, 'lb', uptake, -Inf);
model = setParam(model, 'ub', uptake, Inf);

% minimize substrate uptake
model = setParam(model, 'obj',{'r_2111'}, 0);
model = setParam(model, 'obj', uptake, 1);
sol = solveLP(model,1)

printFluxes(model, sol.x, true);

%%
% Set maximum uptake of carbon source back to 1 mmol/gDCW*hr
model = setParam(model, 'lb', 'r_1714', -1);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model,'obj','r_2111',1);
%model = setParam(model, 'eq', 'r_0959', 0);

% To perform simple FBA, use solveLP
sol = solveLP(model, 1)
% Flux distributions can be displayed by printFluxes:
printFluxes(model, sol.x, true);

%printFluxes(model, sol.x, false);

%% Lipid as objective
% Parameters for simulating batch growth in complex medium
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));

rxnList = findRxnsFromMets(model, {'s_0394','s_0397','s_0458','s_0796',...
 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});

targetsNames = model.rxnNames(getIndexes(model,rxnList,'rxns'));
targetsIDS = model.rxns(getIndexes(model,rxnList,'rxns'));
requiredRxns = targetsIDS(endsWith(targetsNames,'exchange'));

% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);

aminoacidRxns = {'r_1810'; ... % L-glycine
                 'r_1873'; ... % L-alanine
                 'r_1879'; ... % L-arginine
                 'r_1880'; ... % L-asparagine
                 'r_1881'; ... % L-aspartate
                 'r_1883'; ... % L-cysteine
                 'r_1889'; ... % L-glutamate
                 'r_1891'; ... % L-glutamine
                 'r_1893'; ... % L-histidine
                 'r_1897'; ... % L-isoleucine
                 'r_1899'; ... % L-leucine
                 'r_1900'; ... % L-lysine
                 'r_1902'; ... % L-methionine
                 'r_1903'; ... % L-phenylalanine
                 'r_1904'; ... % L-proline
                 'r_1905'; ... % L-serine
                 'r_1911'; ... % L-threonine
                 'r_1912'; ... % L-tryptophan
                 'r_1913'; ... % L-tyrosine
                 'r_1914'};    % L-valine
             
model = setParam(model, 'lb', aminoacidRxns, -0.01);

% set glucose as carbon source and unlimited O2 for aerobic growth
model = setParam(model, 'lb', {'r_1714'}, -2.01);    % glucose
model = setParam(model, 'ub', {'r_1714'}, 0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
%
% set ammonium uptake
model = setParam(model, 'ub', {'r_1654'}, -0.35); %ammonium

% set lipid pseudoreaction as objective
model = setParam(model, 'obj',{'r_2108'}, 1);   
%
sol = solveLP(model,1)
%
printFluxes(model, sol.x, true);
%
%printFluxes(model, sol.x, false);

%%

% Set maximum uptake of carbon source back to 1 mmol/gDCW*hr
model = setParam(model, 'lb', 'r_1714', -1);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model,'obj','r_2108',1);

% To perform simple FBA, use solveLP
sol = solveLP(model, 1)
% Flux distributions can be displayed by printFluxes:
printFluxes(model, sol.x, true);

%printFluxes(model, sol.x, false);

%% 3.9 SAVE TO GITHUB
% The model is tracked and distributed via a GitHub repository.
% Before saving the model, we will add some extra information.
model.annotation.defaultLB    = -1000; % Default lower bound
model.annotation.defaultUB    = +1000; % Default upper bound
model.annotation.taxonomy     = 'taxonomy/460523';
model.annotation.givenName    = 'Rafaela';
model.annotation.familyName   = 'Ventorim';
model.annotation.email        = 'rafaela.ventorim@ufv.br';
model.annotation.organization = 'Universidade Federal de Vicosa';
model.annotation.note         = 'First draft model';
model.id                      = 'papla';
model.description             = 'Papiliotrema laurentii-GEM';

% Save workspace
save([root 'scrap/finalmodel.mat'])
% load([root 'scrap/finalmodel.mat'])

%Extport to inspect:
exportToExcelFormat(model,[root '/scrap/finalmodel_paplaGEM.xlsx']);
exportModel(model,[root 'scrap/finalmodel_paplaGEM.xml'])

%newCommit(model);

