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

%% %% 3 Lipid curation

% Curate lipid metabolism, by modifying reactions with metabolites that
% have an acyl-chain and adding 16:1 lipid chain.
load([root '/scrap/model_r2.mat']);
load([root '/scrap/modelTemplate.mat']); % rhto-GEM 

%% Curate reactions specified in lipidTemplates.txt
% Load text file that contains template reactions, in which compartments
% they are located, which grRules should be applied (specified for each
% compartment), and which sets of acyl-chains are used.

% Acyl chains of biotechnological importance 16:0, 16:1, 18:0, 18:1, 18:2, 
% 18.3. Assume sn1 position to be satured (16:0 and 18:0). Assume sn2 
% position to be unsatured (16:1, 18:1, 18:2, 18:3).

fid         = fopen([data '/reconstruction/lipidTemplates.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', '\t');
fclose(fid);

% Reorganize the content so that it can be used by the addLipidReactions
% function.
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.grRules     = loadedData{4};
template.chains = {};
for k = 1:numCols-3; template.chains(:,k) = loadedData{k+4}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');
% Remove reactions that match a lipid template reaction (ignoring acyl-chains)
toRemove    = regexprep(template.rxns,'CHAIN.*','');
toRemove    = find(startsWith(model.rxnNames,toRemove));
model       = removeReactions(model,toRemove);

% Now use the templates to add the relevant reactions to the model. If a
% reaction already existed in the R. toruloides template model, then it
% will use the same reaction identifier.
model = addLipidReactions(template,model,modelRhto);

%% Add lipid transport reactions
fid         = fopen([data '/reconstruction/lipidTransport.txt']);
firstLine   = fgets(fid);
numCols     = numel(strfind(firstLine,char(9))); % number of \t
loadedData  = textscan(fid, [repmat('%q ', [1, numCols]) '%q'], 'delimiter', '\t');
fclose(fid);

clear template
template.rxns   = loadedData{1};   template.eqns        = loadedData{2};
template.comps  = loadedData{3};   template.chains = {};
for k = 1:numCols-2; template.chains(:,k) = loadedData{k+3}; end
template.chains = regexprep(template.chains,':0(\d)', ':$1');
toRemove    = regexprep(template.rxns,'CHAIN.*','');
toRemove    = find(startsWith(model.rxnNames,toRemove));
model       = removeReactions(model,toRemove);

model = addLipidReactions(template,model,modelRhto);

%% Add SLIME reactions
clear template
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

model = addMetabolite(model, 's_3741','C16:1 chain', 'C16H30O2');
model.metComps(2326) = 1;

model=addSLIMEreactions(template,model,modelRhto);

chainExIdx  = getIndexes(model,'r_4064','rxns');
backbExIdx  = getIndexes(model,'r_4062','rxns');
model       = setParam(model,'ub',[chainExIdx,backbExIdx],1000);
solveLP(model,1)

model           = addTransport(model,'c','erm',{'palmitate','stearate','oleate','linoleate','linolenate','palmitoleate'},true,false,'t_');
model           = addTransport(model,'c','ce',{'palmitate','stearate','oleate','linoleate','linolenate','palmitoleate'},true,false,'t_');

save([root '/scrap/model_r3.mat'],'model');

% To inspect the second draft model:
exportToExcelFormat(model,[root '/scrap/model_r3.xlsx']);

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])
%cd('..'); newCommit(model);cd('reconstruction')

%% %% 4
% DEFINE BIOMASS COMPOSITION

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

save([root 'scrap/model_r4.mat'])
% load([root 'scrap/model_r4.mat'])
clear indexes equations loadedData fid BM biomassRxns ans

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Extport to inspect:
exportToExcelFormat(model, [root 'scrap/model_r4.xlsx']);

%% %% 5 GAP-FILLING
% 5A MENECO
% MENECO requires the target compounds to already be part of the draft
% model. This should be fine here, as we added the whole biomass equation
% above.
load([root '/scrap/model_r4.mat'],'model');
load([root '/scrap/modelTemplate.mat']);

%% Export model for meneco
% exportModel(model,'../meneco/r4_preMeneco.xml')
% exportModel(modelRhto,'../meneco/rhtoRepair.xml')
% exported externally using cobrapy

% Find targets: any substrate for the pseudoreactions, us the following
% text to reconstruct menecoTargets.sbml.
% rxnIdx  = find(contains(model.rxnNames,'pseudoreaction'));
% targets = find(any(model.S(:,rxnIdx)<0,2));
% [model.mets(targets), model.metNames(targets)]
% targetSBML=strcat('<species id="M_',model.mets(targets),...
%     '" name="',model.metNames(targets),'"/>');

% Identified by MENECO (see meneco.txt for output file).
fid         = fopen([data '/meneco/menecoRxns.txt']);
menecoRxns  = textscan(fid,'%s'); fclose(fid);
menecoRxns  = menecoRxns{1};

% If these reactions are present, that means that their respective enzymes
% are present. Any other reaction annotated to the same enzymes should also
% be added.
menecoRxns  = getAllRxnsFromGenes(modelRhto,menecoRxns);
model       = addRxnsGenesMets(model,modelRhto,menecoRxns,true,'Identified by MENECO to produce biomass components',1);

% Test meneco results, using the same exchange reactions as in menecoSeeds.sbml.
% model_tmp=addExchangeRxns(model,'in',{'s_0394','s_0397','s_0458','s_0796',...
% 's_0805','s_1198','s_0420','s_1277','s_1324','s_1468','s_0565','s_0452',...
% 's_0529','s_0925','s_1616','s_1582','s_1583','s_1585','s_1587','s_1589',...
% 's_1590','s_1591','s_1593','s_1594','s_1596','s_1598','s_1600','s_1602',...
% 's_1604','s_1606','s_1607','s_1608','s_1610','s_1612','s_1614'});
% out=cell(length(targets),2);
% for i=1:length(targets)
%     tmpmodel=addExchangeRxns(model_tmp,'out',targets(i));
%     tmpmodel=setParam(tmpmodel,'obj',numel(tmpmodel.rxns),1);
%     out(i,1)=tmpmodel.metNames(getIndexes(tmpmodel,targets(i),'mets'));
%     soltmp=solveLP(tmpmodel,1);
%     out(i,2)=num2cell(soltmp.f);
% end
% out

model   = setParam(model,'obj','r_2111',1);
sol     = solveLP(model,1)

%% 5B Gap-filling RAVEN
% After MENECO, model still couldn't produce biomass. We decided to further
% curate model using RAVEN gap-filling function.

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

% From the Yli model, remove all exchange reactions (the
% necessary ones we already added, don't want to add new ones)
modelYli2 = removeReactions(modelYli, getExchangeRxns(modelYli, 'both'), true, true, true);

% Run fillGaps function
[~, ~, addedRxns, model] = fillGaps(model, {modelRhto2, modelYli2}, false, true);

% Verify that model can now grow
sol = solveLP(model, 1)
printFluxes(model, sol.x, true)

% We find that the model can now produce biomass after RAVEN gap-filling
% function added the following reactions: r_0132, r_0438 and r_4046

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

% Water and H+ exchange was unbalanced in rich media simulations due to
% amino acid reactions from iYali that are not required for biomass 
% production which resulted in dead end fluxes, requiring their removal
model = removeReactions(model,{'y300065','y300066','y200008'});

save([root '/scrap/model_r5.mat'],'model');
% load([root 'scrap/model_r5.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/model_r5.xlsx']);

clear addedRxns rxns sol biomassRxns fid menecoRxns modelRhto2 modelYli2

%% 6 Add and adjust lista-specific reactions

% Manually add a list of metabolites
metsToAdd.metNames      = {'L-Rhamnose','L-Rhamnose','L-rhamnonic acid-gamma-lactone','L-rhamnonic acid','3,6-dideoxy-L-erythro-hexulosonic acid','lactose','lactose','cellobiose','levoglucosan','levoglucosan'};
metsToAdd.compartments  = {'e','c','c','c','c','e','c','e','e','c'};
metsToAdd.metFormulas   = {'C6H12O5','C6H12O5','C6H10O5','C6H12O6','C6H9O5','C12H22O11','C12H22O11','C12H22O11','C6H10O5','C6H10O5'};
metsToAdd.mets          = generateNewIds(model,'mets','s_',length(metsToAdd.metNames));
%metsToAdd.metMiriams    = repmat({struct('name',{{'sbo'}},'value',{{'SBO:0000649'}})},1,3);
model                   = addMets(model,metsToAdd); clear metsToAdd;

% Manually add a list of reactions

fid         = fopen([data '/reconstruction/listaSpecificRxns.txt']);
loadedData  = textscan(fid,'%q %q %q %q %q','delimiter',{'\t','"'},'MultipleDelimsAsOne',1,'TreatAsEmpty','_');
fclose(fid);

clear rxnsToAdd
rxnsToAdd.equations     = regexprep(loadedData{1},'***','');
rxnsToAdd.rxnNames      = regexprep(loadedData{2},'***','');
rxnsToAdd.grRules       = regexprep(loadedData{3},'***','');
rxnsToAdd.subSystems    = regexprep(loadedData{4},'***','');
for i=1:numel(rxnsToAdd.subSystems)
    rxnsToAdd.subSystems{i}=rxnsToAdd.subSystems(i);
end
rxnsToAdd.eccodes       = regexprep(loadedData{5},'***','');
rxnsToAdd.rxns          = generateNewIds(model,'rxns','r_',length(rxnsToAdd.rxnNames));
model                   = addRxns(model,rxnsToAdd,3,'',false,true); clear rxnsToAdd
model                   = removeReactions(model,'r_0719',false,false,false);

save([root '/scrap/model_r6.mat'],'model');
% load([root 'scrap/model_r6.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/model_r6.xlsx']);

clear ans fid i loadedData

%% 7 Update gene-reaction associations

% All Rhto genes have a RHTO in the name, while Lista genes do not.
rxnIdx      = strfind(model.grRules,'RHTO');
rxnIdx      = ~cellfun('isempty',rxnIdx); % Get reaction indexes
out         = cell(length(find(rxnIdx)),3);
out(:,1)    = model.rxns(rxnIdx);
out(:,2)    = model.rxnNames(rxnIdx);
out(:,3)    = model.grRules(rxnIdx);

% From this list, through manual curation define the following new grRules
fid         = fopen([data '/reconstruction/updateGrRules.txt']);
loadedData  = textscan(fid,'%s %s','delimiter','\t'); fclose(fid);
rxns        = loadedData{1};
grRules     = loadedData{2};

model = changeGrRules(model,rxns,grRules,true);

% Correct faulty grRules where the same complex is representated multiple
% times
for n = 1:length(model.grRules)
    if any(model.grRules{n})
        noAnd = strfind(model.grRules(n),'and');
        noAnd = any(vertcat(noAnd{:})); % Give 0 if no 'and' is present.
        if noAnd == 0
            geneList = transpose(cell(unique(regexp(model.grRules{n},'[)(]*|( and )*|( or )*','split'))));
            geneList = regexprep(geneList,'[(*)*]','');
            if length(geneList) == 1
                newgrRule = geneList;
            else
                newgrRule = geneList{1};
                for k = 2:length(geneList)
                    newgrRule = [newgrRule ' or ' geneList{k}];
                end
            end
            model.grRules(n) = cellstr(newgrRule);
        end
    end
end

model = deleteUnusedGenes(model);

save([root '/scrap/model_r7.mat'],'model');
% load([root 'scrap/model_r7.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/model_r7.xlsx']);

clear ans fid grRules loadedData rxns

%% Remove Double NGAM

%load([root 'scrap/model_r7.mat'])

model = removeReactions(model,'r_4046_rhto')

% Remove 'rhto' from subsystems
model.subSystems = cellfun(@(x) regexprep(x,'rhto[0-9]+ +',''),model.subSystems, 'UniformOutput', 0);

% Remove unused metabolites
model = removeMets(model,all(model.S == 0,2),false,true,true,true);

save([root '/scrap/model_r8.mat'],'model');
% load([root 'scrap/model_r8.mat'])

disp(['Number of genes / rxns / mets in model:  ' ...
    num2str(length(model.genes)) ' / ' ...
    num2str(length(model.rxns)) ' / ' ...
    num2str(length(model.mets))])

% Export to Excel format for easy inspection
exportToExcelFormat(model,[root '/scrap/model_r8.xlsx']);


%% Export model:



%Add model information.
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



%Save model
exportForGit(model,'lista-GEM','..',{'mat', 'txt', 'xlsx', 'xml'});

exportModel(model,'lista-GEM.xml')

<<<<<<< HEAD



=======
>>>>>>> main
