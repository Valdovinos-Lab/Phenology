%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for running specialization metrics on alpha values from 
% model output. Specialization metrics include Deviation from
% Generalization (DFG) and Coefficient of Variation (CV). 
% Original Authors: Fernanda S. Valdovinos & Paul Glaum
% Edited by Paul Glaum 
% Last Modification: July, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Pollinator species should be listed in columns and each column 
% should sum 1
% Input: M stands for matrix of alphasf or avgAlphasf
% Outputs: 
% v_dfg = the DFG score for each pollinator 
% v_dfg_st the DFG score for each pollinator standardized to the number of 
%   plant species, range [0,1]. Used in the publication
% v_dfg_Max = the maximum pairwise difference between a single pollinator’s
%   alpha values on different plants
% v_SD = the standard deviation of the pairwise differences in alpha values
% v_CV = the coefficient of variance in the alpha values. Used in the 
%   publication
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_dfg, v_dfg_st, v_dfg_Max, v_SD, v_CV]=DFG(M) 

[rows, cols]=size(M);

for i=1:cols
    Mdiff=abs(minus(M(:,i),M(:,i)'));
    v_dfg(i,1)=sum( sum( triu(Mdiff) ) );
    v_dfg_Max(i,1)=max( max( triu(Mdiff) ) );
    v_SD(i,1)=std(abs(M(:,i)));
    v_CV(i,1)=std(abs(M(:,i)))/mean(abs(M(:,i)));
    
end

v_dfg_st=v_dfg./(rows-1);% DFG standarized by the max deviation, i.e. specialist
                         % with only one link have DFG = number of potential
                         % links - 1.