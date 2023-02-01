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
model_tmp = setParam(model_tmp, 'eq', {'r_1714'}, 0);    % glucose
%model_tmp = setParam(model_tmp, 'lb', {'r_1714'}, -3);
%model_tmp = setParam(model_tmp, 'ub', {'r_1714'}, 0);    

%model_tmp = setParam(model_tmp, 'lb', {'r_1634'}, -1);  % acetate
%model_tmp = setParam(model_tmp, 'lb', {'r_1878'}, -3);  % arabinose
%model_tmp = setParam(model_tmp, 'lb', {'r_4348'}, -3);  % cellobiose
%model_tmp = setParam(model_tmp, 'lb', {'r_1687'}, -3);  % citrate
%model_tmp = setParam(model_tmp, 'lb', {'r_1761'}, -3);  % ethanol
%model_tmp = setParam(model_tmp, 'lb', {'r_1710'}, -3);  % galactose
%model_tmp = setParam(model_tmp, 'lb', {'r_4345'}, -3);  % lactose
model_tmp = setParam(model_tmp, 'lb', {'r_4350'}, -3);  % levoglucosan
%model_tmp = setParam(model_tmp, 'lb', {'r_4339'}, -10);  % rhamnose
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

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\Lipomyces\GEM-Lipomyces\ModelFiles\mat\lista-GEM.mat')

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

ds = 0.05
do = 0.25
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

ds = 0.05
do = 0.25
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

ds = 0.025
do = 0.125
s = 0:ds:5
o = 0:do:25


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

load('C:\Users\dudul\OneDrive\Documentos_UFV\LABFIS\Lipomyces\GEM-Lipomyces\ModelFiles\mat\lista-GEM.mat')

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

ds = 0.05
do = 0.25
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

ds = 0.05
do = 0.25
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

ds = 0.025
do = 0.125
s = 0:ds:5
o = 0:do:25

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

