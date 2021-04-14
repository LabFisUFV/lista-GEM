function newCommit(model)

%Check if in master:
currentBranch = git('rev-parse --abbrev-ref HEAD');
if strcmp(currentBranch,'master')
    error('ERROR: in master branch. For new releases, use newRelease.m')
end

if ~exist('model')
    %Load model:
    model = importModel('../ModelFiles/xml/lista-GEM.xml');
end

%Save model
exportForGit(model,'lista-GEM','../',{'txt', 'xml', 'yml'});
end