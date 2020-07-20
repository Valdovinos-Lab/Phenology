%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for integrating phenological modifiers into plant and
% pollinator behavior and interactions in plant-pollinator networks
% using Valdovinos et al. 2013, 2016 model as basis
% Original Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Modified and Edited by Paul Glaum 
% Last Modification: July, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functions runs a single simulation of the phenology model. 
% Inputs: 
% In: Input network 
% bloomSpan: Time step duration of flowering bloom period
% breakValue: The number of cycles between subsequent bloom periods
% flightSpan: Time step duration of pollinator flight period
% breakValueF: The number of cycels between subsequet flight periods
% Outputs: 
% Densities: p (plants), nectar(rewards), a(animals), alphas(effort)
% Varients: i=initial, f=final, avg...f=avg over the last final 1000 steps
% indiePlantOverlap: Average Resource Overlap (ARO)
% overallPlantOverlap: Total Resource Overlap (TRO)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pI, nectari, ai, alphasi, pf, nectarf, af, alphasf, nectar, avgAlphasf, indiePlantOverlap, overallPlantOverlap]=run_Phenology_model(In,bloomSpan,breakValue,flightSpan,breakValueF)

% Legacy Inputs:
%    ifrG = Fraction of pollinators with adaptive foraging (AF): 0=0, 1=0.25,
%                   2=0.5, 3=0.75, 4=1.
%    muAP = indicates level of mortality of plants and animals without AF or rem:
%    1= animals die, 2=plants die, 3=both live, 4=both die 
%   **Note, results of were produced using both set at 4. 
ifrG=4; muAP=4;

global J_pattern

%Optional seed set for random number generator. 
%rand('seed',2);

[rows cols]= size(In);
J_pattern = J_zero_pattern(In) ;
frG=ifrG*0.25;

vectG=Gprim(frG,cols);
[pI, nectari, ai, alphasi, pf, nectarf, af, alphasf, nectar, avgAlphasf]=Val_ext_phenology(vectG,In,muAP,bloomSpan,breakValue,flightSpan,breakValueF);
pf=full(pf);
nectarf=full(nectarf);
af=full(af);
alphasf=full(alphasf);
nectar=full(nectar);
avgAlphasf=full(avgAlphasf);

m=rows; %the number of plants
n=cols; %the number of pollinators

indiePlantOverlap=zeros(m,1); %storage for resource overlap of individual plants
overallPlantOverlap=0; %overall resource overlap

%Calculate the temporal overlap in resource availability using the process below
for i=1:(m-1)
    
    for j=(i+1):m
        arr=cat(3,nectar(:,i),nectar(:,j) ); %add each Rewards vector to pairwise concatenated array
        overlapData=min(arr,[],3); %create overlapping flowering time between these two flowers
        thisOverlap=trapz(overlapData); %take intregral to get overlap of flowering time
        
        %Add this to quantify total resource overlap
        overallPlantOverlap=overallPlantOverlap+thisOverlap; 
        %Add this to quantify resource overlap of flower i
        indiePlantOverlap(i)=indiePlantOverlap(i)+thisOverlap;
        %Add this to quantify resource overlap of flower j
        indiePlantOverlap(j)=indiePlantOverlap(j)+thisOverlap;
        
    end
end

indiePlantOverlap=full(indiePlantOverlap);
overallPlantOverlap=full(overallPlantOverlap);


end