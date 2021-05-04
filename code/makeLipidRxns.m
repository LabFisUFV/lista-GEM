function [newEqns, newNames] = makeLipidRxns(templEqn, templName, chains, comps)
% makeLipidRxns
%   Given a template reaction equation and name, new equations and names
%   generated using the specified fatty acid chains on the indicated
%   positions, while modifying the compartment.
%
%   templEqn        template reaction equation (see example below)
%   templName       template reaction name (see example below)
%   chains          
%   comps           string or cell array of compartment abbreviations, such
%                   as 'ce' or {'mm','erm'}
%
%   newEqns         cell array of strings with new reaction equations
%   newNames        cell array of strings with new reaction names
%
%   Example of use: assume the following parameters
%   templEqn    = 'CHAIN1_COA[COMP] + gly-3-P[COMP] => CoA[erm] +
%                  1-acyl-gly-3-P (CHAIN1_NUMB)[COMP]'
%   templName   = 'gly-3-p acyltransferase (CHAIN1), COMP'
%   comps       = {'erm','lp'}
%   chain1      = {'16:0','18:3'}
%
%   [newEqns newNames] =
%   makeLipidRxns(templEqn,templName,false,comps,chain)
%
%   will generate four reaction names and equations, two for each
%   COMP compartment and two for each CHAIN fatty acid chain. CHAIN1 to
%   CHAIN4 refer to the fatty acid chain positions (were relevant, in the
%   example there is only one fatty acid chain). _COA specifies that the
%   fatty-acyl-CoA should be used (e.g. palmitoyl-CoA), while _NUMB
%   specifies that the fatty acid should be represented in number style
%   (e.g. 16:0). A third alternative (not shown here), is _ATE, where the
%   fatty acid is presented as e.g. palmitate.
%
%   If 'reduced' is set to true, a reduced set of reactions are given,
%   filtered for duplicate position-aspecific fatty acids. Example:
%   Reaction 1: chain 1 = 16:0, chain 2: 16:1, chain 3: 18:0
%   Reaction 2: chain 1 = 16:1, chain 2: 18:0, chain 3: 16:0
%   Reaction 3: chain 1 = 16:0, chain 2: 18:0, chain 3: 16:0
%   Both reaction 1 and 2 involve 1x 16:0, 1x 16:1 and 1x 18:0. Reaction 2
%   will therefore be filtered out. This reduces the number of possible
%   combinations, but retains the overall profile.
%
%   Compartment and fatty acid names and abbreviations might have to be
%   adjusted for your model, this is done in the code below, under fillEqns
%   and fillComps.
%
%   Eduard Kerkhoven, 2019-06-25


chainDB.numb     = {'16:0', '16:1', '18:0', '18:1', '18:2', '18:3'};
chainDB.ate      = {'palmitate', 'palmitoleate', 'stearate', 'oleate', ...
    'linoleate', 'linolenate'};
chainDB.coa      = {'palmitoyl-CoA', 'palmitoleoyl-CoA', ...
    'stearoyl-CoA', 'oleoyl-CoA', 'linoleoyl-CoA', 'linolenoyl-CoA'};

chains(cellfun('isempty',chains)) = [];

%% Duplicate number of reactions for number of chain compositions
newEqns = cellstr(repmat(templEqn, length(chains),1));
newNames = cellstr(repmat(templName, length(chains),1));

for i = 1:numel(newEqns)
    chainComb = strtrim(split(chains{i},','));
    for j = 1:numel(chainComb)
        POS = ['CHAIN' num2str(j)];
        [~, chainLoc]   = ismember(chainComb{j},chainDB.numb);
        chainAtes       = chainDB.ate(chainLoc);
        chainCoas       = chainDB.coa(chainLoc);
        
        newEqns(i) = regexprep(newEqns(i),[POS '_NUMB'],chainComb{j});
        newEqns(i) = regexprep(newEqns(i),[POS '_ATE'],chainAtes);
        newEqns(i) = regexprep(newEqns(i),[POS '_COA'],chainCoas);
        newNames(i)= regexprep(newNames(i),[POS '_COA'],chainCoas);
        newNames(i)= regexprep(newNames(i),[POS '_ATE'],chainAtes);
        newNames(i)= regexprep(newNames(i),POS,chainComb{j});
    end
end

if exist('comps')
    [newEqns, newNames] = fillComps(newEqns, newNames, comps);
end

end

%% Duplicate for compartments
function [newEqns, newNames] = fillComps(templEqns, templNames, comps)

compsDB.id      = {'c','ce','e','er','erm','g','gm','lp','m','mm','n',...
    'p','v','vm'};
compsDB.name      = {'cytoplasm','cell envelope','extracellular',...
    'ER reticulum','ER membrane','Golgi',...
    'Golgi membrane','lipid particle','mitochondrion',...
    'mitochondrial membrane','nucleus','peroxisome','vacuole',...
    'vacuolar membrane'};

if ~iscell(templEqns)
    templEqns   = {templEqns};
end

[~, compLoc]    = ismember(comps,compsDB.id);
compNames       = compsDB.name(compLoc);

if ~iscell(comps)
    comps       = {comps};
    newEqns     = cell(length(templEqns),1);
    newEqns(:)  = templEqns;
    newNames    = cell(length(templNames),1);
    newNames(:) = templNames;    
else
    newEqns     = repmat(templEqns,length(comps),1);
    newNames    = repmat(templNames,length(comps),1);    
end

idx = 1:1:length(newEqns);
idx = reshape(idx,[length(idx)/length(comps),length(comps)]);

for i=1:length(comps)
    newEqns(idx(:,i))  = regexprep(newEqns(idx(:,i)),'\[COMP\]',['\[' comps{i} '\]']);
    newNames(idx(:,i)) = regexprep(newNames(idx(:,i)),'COMP',['' compNames{i} '']);
end
end
