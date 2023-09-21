%%SIMULATIONS

clear model_tmp
model_tmp = model;

% Parameters for simulating batch growth in minimal medium
exchangeRxns = model_tmp.rxns(endsWith(model_tmp.rxnNames,'exchange'));

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange

% block all uptake and allow only required metabolites
model_tmp = setParam(model_tmp, 'lb', exchangeRxns, 0);
model_tmp = setParam(model_tmp, 'lb', requiredRxns, -1000);
%
% set glucose as carbon source and unlimited O2 for aerobic growth
%model_tmp = setParam(model_tmp, 'eq', {'r_1714'}, 0);    % glucose
model_tmp = setParam(model_tmp, 'lb', {'r_1714'}, -3);
%model_tmp = setParam(model_tmp, 'ub', {'r_1714'}, 0);    

%model_tmp = setParam(model_tmp, 'lb', {'r_1634'}, -3);  % acetate
%model_tmp = setParam(model_tmp, 'lb', {'r_1878'}, -3);  % arabinose
%model_tmp = setParam(model_tmp, 'lb', {'r_4348'}, -3);  % cellobiose
%model_tmp = setParam(model_tmp, 'lb', {'r_1687'}, -3);  % citrate
%model_tmp = setParam(model_tmp, 'lb', {'r_1761'}, -3);  % ethanol
%model_tmp = setParam(model_tmp, 'lb', {'r_1710'}, -3);  % galactose
%model_tmp = setParam(model_tmp, 'lb', {'r_4345'}, -3);  % lactose
%model_tmp = setParam(model_tmp, 'lb', {'r_4350'}, -3);  % levoglucosan
%model_tmp = setParam(model_tmp, 'lb', {'r_4339'}, -3);  % rhamnose
%model_tmp = setParam(model_tmp, 'lb', {'r_1546'}, -3);  % R-lactate
%model_tmp = setParam(model_tmp, 'lb', {'r_1551'}, -3);  % S-lactate
%model_tmp = setParam(model_tmp, 'lb', {'r_1715'}, -3);  % mannose
%model_tmp = setParam(model_tmp, 'lb', {'r_1650'}, -3);  % trehalose
%model_tmp = setParam(model_tmp, 'lb', {'r_2104'}, -3);  % xylitol

% restriction of other reactions
model_tmp = setParam(model_tmp, 'lb', {'r_1992'}, -1000);    % O2
model_tmp = setParam(model_tmp, 'ub', {'r_1992'}, 0);

%model_tmp = setParam(model_tmp, 'eq', {'r_1671'}, 0);    % biotin

% set biomass pseudoreaction as objective
model_tmp = setParam(model_tmp, 'lb', {'r_2111'}, 0);   % block biomass uptake
model_tmp = setParam(model_tmp, 'obj',{'r_2111'}, 1);   

%
sol = solveLP(model_tmp, 1)
%
printFluxes(model_tmp, sol.x, true);
%
%printFluxes(model_tmp, sol.x, false);

%%
% nutrient uptake reactions to simulate rich medium conditions
clear model_tmpYPD
model_tmpYPD = model_tmp;

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
              
model_tmpYPD = setParam(model_tmpYPD, 'lb', aminoacidRxns, -1);

solYPD = solveLP(model_tmpYPD)

printFluxes(model_tmpYPD, solYPD.x, true);
% printFluxes(model_tmpYPD, solYPD.x, false);

%% Chemostat simulation
% adjust parameters for chemostat growth
clear
load('matlab_r5.mat')

model = setParam(model, 'eq', {'r_2111'}, 0.06);  % fix specific growth rate at the dilution rate value

uptake = find(strcmp(model.rxnNames,'D-xylose exchange')); % remove constraints on substrate uptake
model = setParam(model, 'lb', uptake, -Inf);
model = setParam(model, 'ub', uptake, Inf);

model = setParam(model, 'eq', {'r_1714'}, 0);

model = setParam(model, 'lb', {'r_2091'}, -Inf);
model = setParam(model, 'eq', {'r_1654'}, 0);

% minimize substrate uptake
model = setParam(model, 'obj',{'r_2111'}, 0);
model = setParam(model, 'obj', uptake, 1);
sol = solveLP(model, 1)

printFluxes(model, sol.x, true);

%%
clear model_CC
model_CC = model_tmp;

%model_CC = setParam(model_CC, 'lb', {'r_2091'}, -0.025);
%model_CC = setParam(model_CC, 'eq', {'r_1654'}, 0);

model_CC = setParam(model_CC, 'eq', {'r_1714'}, 0);

solution = simulateChemostat(model_CC, 1, ['r_1718' 'r_2111']);
printFluxes(model_CC, solution, true);

%% Surface substrate vs. oxygen for growth
%Minimal media 

load('../model/mat/lista-GEM.mat')

model.annotation.defaultLB    = -1000; % Default lower bound
model.annotation.defaultUB    = +1000; % Default upper bound

model_tmp = model;

% Parameters for simulating batch growth in minimal medium
exchangeRxns = model_tmp.rxns(endsWith(model_tmp.rxnNames,'exchange'));

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange

% block all uptake and allow only required metabolites
model_tmp = setParam(model_tmp, 'lb', exchangeRxns, 0);
model_tmp = setParam(model_tmp, 'lb', requiredRxns, -1000);
model_tmp = setParam(model_tmp, 'lb', {'r_2111'}, 0);   % block biomass uptake
model_tmp = setParam(model_tmp, 'obj', {'r_2111'}, 1);   % biomass
%%
%%
%glucose 

ds = 0.20
do = 1.0
s = 0:ds:10
o = 0:do:50


for i = 1:length(s)
    for j = 1:length(o)
model_tmp = setParam(model_tmp,'lb', 'r_1714',-s(i));
model_tmp = setParam(model_tmp,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model_tmp,'max');
growthRates(i,j) = FBAsolution.f;
    end
end

pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot


%xylose

model_tmp = setParam(model_tmp,'eq', 'r_1714',0);

ds = 0.20
do = 1.0
s = 0:ds:10
o = 0:do:50



for i = 1:length(s)
    for j = 1:length(o)
model_tmp = setParam(model_tmp,'lb', 'r_1718',-s(i));
model_tmp = setParam(model_tmp,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model_tmp,'max');
growthRates(i,j) = FBAsolution.f;
    end
end


pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot


%glycerol
model_tmp = setParam(model_tmp,'eq', 'r_1714',0);

ds = 0.20
do = 1.0
s = 0:ds:10
o = 0:do:50


for i = 1:length(s)
    for j = 1:length(o)
model_tmp = setParam(model_tmp,'lb', 'r_1808',-s(i));
model_tmp = setParam(model_tmp,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model_tmp,'max');
growthRates(i,j) = FBAsolution.f;
    end
end



pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot




%% Surface substrate vs.oxygen for TAG production
% Minimal media 

load('../model/mat/lista-GEM.mat')

idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]', ...
    'triglyceride (1-16:0, 2-18:1, 3-18:1)[lp]'}, 'metcomps');
% Add exchange reactions for products
rxnsToAdd.rxns          = 'exch_TAG';
rxnsToAdd.mets          = model.mets(idx);
rxnsToAdd.stoichCoeffs  = {[-1, -1]}; 
rxnsToAdd.lb            = 0;
model = addRxns(model,rxnsToAdd);

model = setParam(model,'obj','exch_TAG',1);
model = setParam(model,'obj','r_2111',0);
model = setParam(model,'lb','r_4046',0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
model = setParam(model,'eq','r_2111',0.01);
model = setParam(model,'obj','exch_TAG',1);

%glucose 


ds = 0.20
do = 1.0
s = 0:ds:10
o = 0:do:50


for i = 1:length(s)
    for j = 1:length(o)
model = setParam(model,'lb', 'r_1714',-s(i));
model = setParam(model,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model,'max');
growthRates(i,j) = FBAsolution.f;
    end
end


printFluxes(model, FBAsolution.x)

pcolor(o,s,growthRates) %2D plot

%xylose

model = setParam(model,'eq', 'r_1714',0);

ds = 0.20
do = 1.0
s = 0:ds:10
o = 0:do:50


for i = 1:length(s)
    for j = 1:length(o)
model = setParam(model,'lb', 'r_1718',-s(i));
model = setParam(model,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model,'max');
growthRates(i,j) = FBAsolution.f;
    end
end


pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot


%glycerol
model = setParam(model,'eq', 'r_1714',0);

ds = 0.20
do = 1.0
s = 0:ds:10
o = 0:do:50

for i = 1:length(s)
    for j = 1:length(o)
model = setParam(model,'lb', 'r_1808',-s(i));
model = setParam(model,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model,'max');
growthRates(i,j) = FBAsolution.f;
    end
end



pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot

%% Surface nitrogen vs. oxygen for growth

%Fixado fonte de carbono (-3)
%Minimal media 

load('../model/mat/lista-GEM.mat')

model.annotation.defaultLB    = -1000; % Default lower bound
model.annotation.defaultUB    = +1000; % Default upper bound

model_tmp = model;

% Parameters for simulating batch growth in minimal medium
exchangeRxns = model_tmp.rxns(endsWith(model_tmp.rxnNames,'exchange'));

requiredRxns = {'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange

% block all uptake and allow only required metabolites
model_tmp = setParam(model_tmp, 'lb', exchangeRxns, 0);
model_tmp = setParam(model_tmp, 'lb', requiredRxns, -1000);
model_tmp = setParam(model_tmp, 'lb', {'r_2111'}, 0);   % block biomass uptake
model_tmp = setParam(model_tmp, 'obj', {'r_2111'}, 1);   % biomass
%%

%glucose

model_tmp = setParam(model_tmp,'eq', 'r_1714',-3);

ds = 0.18
do = 0.54
s = 0:ds:9
o = 0:do:27

tic
for i = 1:length(s)
    for j = 1:length(o)
model_tmp = setParam(model_tmp,'lb', 'r_1654',-s(i));
model_tmp = setParam(model_tmp,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model_tmp,'max');
growthRates(i,j) = FBAsolution.f;
    end
end
toc

pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot

%xylose

model_tmp = setParam(model_tmp,'eq', 'r_1718',-3);
model_tmp = setParam(model_tmp,'eq', 'r_1714',0);

ds = 0.18
do = 0.54
s = 0:ds:9
o = 0:do:27

tic
for i = 1:length(s)
    for j = 1:length(o)
model_tmp = setParam(model_tmp,'lb', 'r_1654',-s(i));
model_tmp = setParam(model_tmp,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model_tmp,'max');
growthRates(i,j) = FBAsolution.f;
    end
end
toc


pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot


%glycerol
model_tmp = setParam(model_tmp,'eq', 'r_1714',0);
model_tmp = setParam(model_tmp,'eq', 'r_1808',-3);


ds = 0.18
do = 0.54
s = 0:ds:9
o = 0:do:27

tic
for i = 1:length(s)
    for j = 1:length(o)
model_tmp = setParam(model_tmp,'lb', 'r_1654',-s(i));
model_tmp = setParam(model_tmp,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model_tmp,'max');
growthRates(i,j) = FBAsolution.f;
    end
end
toc

pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot

%% Surface oxygen vs. nitrogen for TAG production
% Minimal media 

load('../model/mat/lista-GEM.mat')

idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]', ...
    'triglyceride (1-16:0, 2-18:1, 3-18:1)[lp]'}, 'metcomps');
% Add exchange reactions for products
rxnsToAdd.rxns          = 'exch_TAG';
rxnsToAdd.mets          = model.mets(idx);
rxnsToAdd.stoichCoeffs  = {[-1, -1]}; 
rxnsToAdd.lb            = 0;
model = addRxns(model,rxnsToAdd);

model = setParam(model,'obj','exch_TAG',1);
model = setParam(model,'lb','r_4046',0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);
model = setParam(model,'eq','r_2111',0.01);

%glucose

model = setParam(model,'eq', 'r_1714',-3);

ds = 0.18
do = 0.54
s = 0:ds:9
o = 0:do:27

tic
for i = 1:length(s)
    for j = 1:length(o)
model = setParam(model,'lb', 'r_1654',-s(i));
model = setParam(model,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model,'max');
growthRates(i,j) = FBAsolution.f;
    end
end
toc

pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot

%xylose

model = setParam(model,'eq', 'r_1718',-3);
model = setParam(model,'eq', 'r_1714', 0);

ds = 0.18
do = 0.54
s = 0:ds:9
o = 0:do:27

tic
for i = 1:length(s)
    for j = 1:length(o)
model = setParam(model,'lb', 'r_1654',-s(i));
model = setParam(model,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model,'max');
growthRates(i,j) = FBAsolution.f;
    end
end
toc

pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot

%glycerol

model = setParam(model,'eq', 'r_1808',-3);
model = setParam(model,'eq', 'r_1714', 0);

ds = 0.18
do = 0.54
s = 0:ds:9
o = 0:do:27


tic
for i = 1:length(s)
    for j = 1:length(o)
model = setParam(model,'lb', 'r_1654',-s(i));
model = setParam(model,'lb', 'r_1992',-o(j));
FBAsolution = optimizeCbModel(model,'max');
growthRates(i,j) = FBAsolution.f;
    end
end
toc

pcolor(o,s,growthRates) %2D plot

surfl(o,s,growthRates) %3D plot



%% eMOMA

%clear;clc;

load('../model/mat/lista-GEM.mat')

model_tmp = model
% Parameters for simulating batch growth in minimal medium
exchangeRxns = model_tmp.rxns(endsWith(model_tmp.rxnNames,'exchange'));

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange

% block all uptake and allow only required metabolites
model_tmp = setParam(model_tmp, 'lb', exchangeRxns, 0);
model_tmp = setParam(model_tmp, 'lb', requiredRxns, -1000);


model_tmp = setParam(model_tmp, 'lb', {'r_1992'}, -1000);    % O2
model_tmp = setParam(model_tmp, 'ub', {'r_1992'}, 0);

model_tmp = setParam(model_tmp, 'lb', {'r_2111'}, 0);   % block biomass uptake
model_tmp = setParam(model_tmp, 'obj',{'r_2111'}, 1);   



% Add TAG exchange reaction
idx = getIndexes(model_tmp, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]'}, 'metcomps');
model_tmp = addExchangeRxns(model_tmp, 'out', idx);

% Fix NGAM to low value
model_tmp = setParam(model_tmp,'eq','r_4046', 0.5);

% Block 
model_tmp = setParam(model_tmp,'eq','r_1761', 0); % ethanol
model_tmp = setParam(model_tmp,'eq','r_1650', 0); % trehalose
model_tmp = setParam(model_tmp,'eq','r_1549', 0); % butanediol 
model_tmp = setParam(model_tmp,'eq','r_2033', 0); % piruvate

% Block exchange reactions in the TCA cycle and others based on literature
% information

model_tmp = setParam(model_tmp, 'eq', {'r_1798', ... % Fumarate
                                       'r_1586', ... % 2-oxoglutarate                                       
                                       'r_1552', ... % malate
                                       'r_1989', ... % oxaloacetate
                                       'r_1815', ... % glyoxylate
                                       'r_1634' ... % acetate
                                       }, 0); 

% Block various other lipids, promote TAG accumulation
model_tmp = setParam(model_tmp, 'eq', {'r_1727', ... % Decanoate
                                       'r_1993', ... % Palmitate
                                       'r_1994', ... % Palmitoleate
                                       'r_2189', ... % Oleate
                                       }, 0); 


%Sterols
model_tmp = setParam(model_tmp, 'eq', {'r_2134', ... % 14-demethyllanosterol
                                       'r_1753', ... % episterol
                                       'r_1757', ... % ergosterol
                                       'r_1788', ... % fecosterol
                                       'r_1915', ... % lanosterol
                                       'r_2106', ... % zymosterol
                                       'r_2137', ... % ergosta-5,7,22,24(28)-tetraen-3beta-ol
                                       }, 0);


model_tmp = setParam(model_tmp, 'obj',{'r_2111'}, 1); % Growth as objective
model_tmp = setParam(model_tmp, 'eq', {'r_1808'}, -3); % glycerol


% Check that model still functions
sol   = solveLP(model_tmp,1) 
printFluxes(model_tmp, sol.x, true);

% Get reference flux distribution by FBA under non-restricted growth conditions

modelRef = model_tmp;
solRef = solveLP(modelRef,1);
printFluxes(modelRef, solRef.x);

% Block nitrogen exchange to mimick nitrogen restricted conditions
modelLim    = setParam(model_tmp, 'eq', 'r_1654', 0);
solLim      = solveLP(modelLim,1); % Confirm that no growth is possible
printFluxes(modelLim, solLim.x);

solMOMA     = MOMA(modelRef,modelLim);
printFluxes(modelRef,solMOMA.x)
solMOMA.x(end)

%% Run eMOMA loop

% Reduce the number of reactions to be tested. Reaction should carry flux
% in at least the FBA of reference condition, and/or the N-restricted MOMA,
% otherwise KO or OE would not have an effect.

nonZeroFlux = find(solMOMA.x ~= 0 | solRef.x ~= 0);

% RefMuKO, LimProdKO, RefMuOE, LimProdOE
out=zeros(numel(nonZeroFlux),5);
f=waitbar(0,'Running eMOMA...');

for i=1:numel(nonZeroFlux)
    waitbar(i/numel(nonZeroFlux),f,sprintf('Running eMOMA... %.1f%%',(i/numel(nonZeroFlux))*100))
    % Knockout
    j=nonZeroFlux(i); % Index of reaction in the model
    modelRefKO  = setParam(modelRef,'eq',j,0);
    try
        solRefKO    = solveLP(modelRefKO,1);
        out(i,1)    = -solRefKO.f; % Growth rate in reference conditions
    catch
        out(i,1)    = 0;
    end
    modelLimKO  = setParam(modelRefKO, 'eq', {'r_1654'}, 0);
    try
        solLimKO    = MOMA(modelRefKO, modelLimKO);
        out(i,2)    = solLimKO.x(end);
    catch
        out(i,2)    = 0;
    end

% Overexpression (forcing 2x higher flux, catch non-functional models)
    noGrowth = true;
    OEfactor = 2.1;
    while (noGrowth==true && OEfactor>1)
        OEfactor = OEfactor - 0.1; % To start at 2, in 10% steps
        modelRefOE  = setParam(modelRef,'eq',j,OEfactor*solRef.x(j));
        try
            solRefOE    = solveLP(modelRefOE,1);
            out(i,3)    = -solRefOE.f; % Growth rate in reference conditions
            noGrowth    = false;
            out(i,5)    = OEfactor;
        catch
            out(i,3)    = 0;
            out(i,5)    = OEfactor;
        end
    end
    modelLimOE  = setParam(modelRefOE, 'eq', {'r_1654'}, 0);
    try
        solLimOE    = MOMA(modelRefOE, modelLimOE);
        out(i,4)    = solLimOE.x(end);
    catch
        out(i,4)    = 0;
    end
end
close(f);


%% Filter to growth rate and TAG exchange from reference condition
refMu = -solRef.f; % Reference growth rate
refEx = solMOMA.x(end); % Reference TAG production in N-restriction

filtRes = out;
filtRes(:,[1,3]) = filtRes(:,[1,3])/refMu; % Normalize growth rates to reference
filtRes(:,[2,4]) = filtRes(:,[2,4])/refEx; % Normalize TAG production to reference

% Only keep those cases where N-restriction resulted in TAG production that
% was at least 5% higher than in the non-mutated strain in N-restriction,
% and the growth rate when N-exchange was allowed was at least 90% of the
% growth rate of the non-mutated strain when allowing N-exchange. Keep
% those reactions where the above is true for either the knockout or the
% overexpression (or both).

prodImpr = 1.02; % Minimum 2% production increase during N-restriction
growRed = 0.90; % Minimum 90% of non-mutated growth rate when N-uptake is allowed
exUp = find((filtRes(:,2)>prodImpr & filtRes(:,1)>growRed) | (filtRes(:,4)>prodImpr & filtRes(:,3)>growRed));
exUp = [nonZeroFlux(exUp), filtRes(exUp,:)];
exUp(:,7) = max(exUp(:,[3,5]),[],2); % Add extra column: maximum production for each reaction
[~,I]=sort(exUp(:,7),1,'descend'); % Sort by maximum production
exUp = exUp(I,:); % Sort by maximum production

exUp = [modelRef.rxns(exUp(:,1)),modelRef.rxnNames(exUp(:,1)),num2cell(exUp)];

% Columns of table now refer to:
% 1: Reaction identifier
% 2: Reaction name
% 3: Position in output vector (not really useful anymore, since we already
%    extracted the reaction identifiers and names
% 4: Growth rate when this reaction is knocked out.
% 5: TAG production when this reaction is knocked out and no N-uptake.
% 6: Growth rate when this reaction is overexpressed.
% 7: TAG production when this reaction is overexpressed and no N-uptake.
% 8: Fraction by which the reaction is overexpressed. Two-fold by default,
%    but could be lower if two-fold did not allow growth.
% 9: Maximum of row 5 and 7, to sort the most promising reactions.

fid = fopen('../data/results/eMOMA_glycerol_noSterolExch.tsv','w');
fprintf(fid, '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n',["rxnID" "rxnName" "idx" ...
    "GR_KO" "EX_KO" "GR_OE" "EX_OE" "OEfactor" "EXmax"]);
for j=1:length(exUp)
    fprintf(fid, '%s\t%s\t%i\t%5.4f\t%5.4f\t%5.4f\t%5.4f\t%i\t%5.4f\n',exUp{j,:});
end
fclose(fid);


%% This script predicts metabolic engineering targets for increased TAG production

%FSEOF is based on the principle that increased production requires a redirection of flux, 
%from originally going toward biomass generation, to ideally going (partially) toward our product of interest.
%The slope parameter derived from FSEOF is indicative of how strong each reaction is contributing toward a shift 
%from growth toward production of the target compound,
%and thereby suggests which reactions are promising targets for overexpression to increase productivity

% Minimal media 

load('../model/mat/lista-GEM.mat')

% Parameters for simulating batch growth in minimal medium
exchangeRxns = model.rxns(endsWith(model.rxnNames,'exchange'));

requiredRxns = {'r_1654'; ... % ammonium exchange
                'r_1832'; ... % H+ exchange
                'r_1861'; ... % iron(2+) exchange
                'r_2005'; ... % phosphate exchange
                'r_2020'; ... % potassium exchange
                'r_2060'}; ... % sulphate exchange

% block all uptake and allow only required metabolites
model = setParam(model, 'lb', exchangeRxns, 0);
model = setParam(model, 'lb', requiredRxns, -1000);

model = setParam(model,'lb','r_4046', 0);
model = setParam(model, 'lb', {'r_1992'}, -1000);    % O2
model = setParam(model, 'ub', {'r_1992'}, 0);

% Add exchange reactions for triglyceride (16:0/18:1/18:1-TAG as target).
idx = getIndexes(model, {'triglyceride (1-16:0, 2-18:1, 3-18:1)[erm]'}, 'metcomps');

% Add exchange reactions for products
model = addExchangeRxns(model, 'out', idx);
% Keep track of ids of exchange reactions
rxn = model.rxns(end);

% Perform FSEOF for TAG on glucose
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [-1, 0, 0]);

targets{1} = FSEOF(model, 'r_2111', rxn, 10, 0.9);

% Perform FSEOF for TAG on xylose
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, -1, 0]);

targets{2} = FSEOF(model, 'r_2111', rxn, 10, 0.9);

% Perform FSEOF for TAG on glycerol
model       = setParam(model, 'eq', {'r_1714', 'r_1718', 'r_1808'}, [0, 0, -1]);

targets{3} = FSEOF(model, 'r_2111', rxn, 10, 0.9);

% Summarize results in table
geneAssoc = ~cellfun('isempty',model.grRules);
for i=1:size(targets,2)
    target(:,i)=targets{i}.logical;
end
for i=1:size(targets,2)
    slope(:,i)=targets{i}.slope;
    slope(~target(:,i),i)=nan;
end

target  = find(sum(target,2) & geneAssoc);
[~,I]=sort(sum(slope(target,:),2,'omitnan'),'descend');
out     = [num2cell(slope(target(I),:)), model.rxnNames(target(I)), model.grRules(target(I))];

fid = fopen('/data/results/fseof_TAG_glu_xyl_gly.csv','w');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\n',["glu_TAG" "xyl_TAG" "gly_TAG" ...
    "rxnName" "grRule"]);
for j=1:length(I)
    fprintf(fid,'%2.2f\t%2.2f\t%2.2f\t%s\t%s\n',out{j,:});
end
fclose(fid);
