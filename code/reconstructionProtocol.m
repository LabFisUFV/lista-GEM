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

%% %% This script generates a first draft of a genome-scale model (here: L.
% starkeyi) using the getModelFromHomology function in RAVEN. It uses
% R. toruloides as a template model.

% BLAST against Rhto genome. L. starkeyi protein fasta obtained from NCBI (código)
blastRhto  = getBlast('lista',[data '..\data\genomes\lista.fasta'],'rhto',[data '..\data\genomes\rhto.fasta']);
blastYli   = getBlast('rhto',[data '..\data\genomes\lista.fasta'],'yli',[data '..\data\genomes\iyali.fasta']);

mkdir([root '/scrap'])

save([root '/scrap/blastStruct.mat'],'blast*');
% load([root '/scrap/blastStruct.mat']);
































