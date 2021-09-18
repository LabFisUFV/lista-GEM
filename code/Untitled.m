%% 3.8 SIMULATIONS

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
%model_tmp = setParam(model_tmp, 'lb', {'r_1714'}, -2.5);
%model_tmp = setParam(model_tmp, 'ub', {'r_1714'}, 0);    

%model_tmp = setParam(model_tmp, 'lb', {'r_1634'}, -2.5);  % acetate
%model_tmp = setParam(model_tmp, 'lb', {'r_1878'}, -2.5);  % arabinose
%model_tmp = setParam(model_tmp, 'lb', {'r_4348'}, -2.5);  % cellobiose
%model_tmp = setParam(model_tmp, 'lb', {'r_1687'}, -2.5);  % citrate
%model_tmp = setParam(model_tmp, 'lb', {'r_1761'}, -2.5);  % ethanol
%model_tmp = setParam(model_tmp, 'lb', {'r_1710'}, -2.5);  % galactose
%model_tmp = setParam(model_tmp, 'lb', {'r_4345'}, -2.5);  % lactose
model_tmp = setParam(model_tmp, 'lb', {'r_4350'}, -2.5);  % levoglucosan
%model_tmp = setParam(model_tmp, 'lb', {'r_4339'}, -2.5);  % rhamnose
%model_tmp = setParam(model_tmp, 'lb', {'r_1546'}, -2.5);  % R-lactate
%model_tmp = setParam(model_tmp, 'lb', {'r_1551'}, -2.5);  % S-lactate
%model_tmp = setParam(model_tmp, 'lb', {'r_1715'}, -2.5);  % mannose
%model_tmp = setParam(model_tmp, 'lb', {'r_1650'}, -2.5);  % trehalose
%model_tmp = setParam(model_tmp, 'lb', {'r_2104'}, -2.5);  % xylitol
%model_tmp = setParam(model_tmp, 'lb', {'r_1808'}, -2.5);  % glycerol
%model_tmp = setParam(model_tmp, 'lb', {'r_1718'}, -2.5);  % xylose
%model_tmp = setParam(model_tmp, 'lb', {'r_2058'}, -2.5);  % sucrose
%model_tmp = setParam(model_tmp, 'lb', {'r_4043'}, -2.5);  % raffinose
%model_tmp = setParam(model_tmp, 'lb', {'r_1931'}, -2.5);  % maltose
%model_tmp = setParam(model_tmp, 'lb', {'r_2056'}, -2.5);  % succinate

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
%clear exchangeRxns requiredRxns model_tmp sol solYPD model_tmpYPD aminoacidRxns
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
              
model_tmpYPD = setParam(model_tmpYPD, 'lb', aminoacidRxns, -0.1);

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