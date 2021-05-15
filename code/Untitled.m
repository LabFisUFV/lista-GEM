%% 3.8 SIMULATIONS

% reducedModel=removeReactions(model,{'y300065','y300066','y200008'});
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
% model_tmp = setParam(model_tmp, 'eq', {'r_1714'}, 0);    % glucose

model_tmp = setParam(model_tmp, 'eq', {'r_1654'}, 0);

model_tmp = setParam(model_tmp, 'lb', {'r_1714'}, -1);
model_tmp = setParam(model_tmp, 'ub', {'r_1714'}, 0);    

model_tmp = setParam(model_tmp, 'lb', {'r_2091'}, -1000);    % urea

model_tmp = setParam(model_tmp, 'lb', {'r_1992'}, -1.890049);    % O2
model_tmp = setParam(model_tmp, 'ub', {'r_1992'}, 0);
%
% set biomass pseudoreaction as objective
model_tmp = setParam(model_tmp, 'lb', {'r_2111'}, 0);   % block biomass uptake
%model_tmp = setParam(model_tmp, 'ub', {'r_2111'}, 1000);
model_tmp = setParam(model_tmp, 'obj',{'r_2111'}, 1);   
%
sol = solveLP(model_tmp, 1)
%
printFluxes(model_tmp, sol.x, true);
%
%printFluxes(model_tmp, sol.x, false);

%%
% nutrient uptake reactions to simulate complex medium conditions
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