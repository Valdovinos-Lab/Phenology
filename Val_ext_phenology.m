%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for integrating phenological modifiers into plant and
% pollinator behavior and interactions in plant-pollinator networks
% using Valdovinos et al. 2013, 2016 model as basis
% Original Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Modified and Edited by Paul Glaum 
% Last Modification: July, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function generates parameters of the dynamic model and 
% safe them in the global structure 'network_metadata' (using
% 'create_metadata.m') so they can be called repeatedly, and integrates the
% ordinary differential equations of the model specified in
% polinizacion_rhs_phenology.m. Once the model is run, results are ploted
% with the time series of plants (top), rewards (middle) and animals 
% (bottom). Each line in the figure represent the density of a species 
% (y-axis) over time (x-axis).
% Inputs: 
% vectG: Legacy input defines which pollinator species exhibit adaptive foraging.
% muAP: Legacy input defining the mortality scenario
% Other inputs are transfered from run_Phenology_model.m 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pI, nectari, ai, alphasi, pf, nectarf, af, alphasf, nectar, avgAlphasf]=Val_ext_phenology(vectG,In,muAP,bloomSpan,breakValue,flightSpan,breakValueF)

global network_metadata J_pattern %indRemP %indRemA %

[m, n]=size(In);
B=sparse(In) ;

% Parameters for uniform distributions
varp=1e-3;%1e-1, 1e-2 %last changed from ...
vara=5e-3;%1e-4, 1e-1 %last changed from 1e-1
%conversion efficiency 
mC=0.2; vC=vara;
%kapa? half saturation constant
mk=0.6; vk=vara;
%seeds per pollination event
mE=0.8; vE=varp;
%pollinator extraction efficiency
mb=0.4; vb=vara;
%interspecific comp coefficient
mU=0.002; vU=varp;
%intraspecific comp coefficient
mw=1.2; vw=varp;
%beta, per individual resources production rate
mB=0.2; vB=varp;
%adaptive rate of foraging, G=2
mG=2; vG=vara;
%background recruitment from seeds to plants, mg=.4
mg=0.4; vg=varp;
%floral resources self limitation parameter, mphi=.04
mphi=0.04; vphi=varp;
vmA=vara;
vmP=varp;

if muAP==1
    mmA=0.01; mmP=0.002; % % high pollinator mortality
elseif muAP==2
    mmA=0.001; mmP=0.02; % % high plant mortality
elseif muAP==3
    mmA=0.001; mmP=0.002; % low plant and animal mortality
else
    mmA=0.003; mmP=0.001; % mortalities used in this published work 
end

%Values for the model parameters obtained from a unifrom distribution
c=uniform_rand(mC,vC,m,n).*B;
k=uniform_rand(mk,vk,m,n).*B;
e=uniform_rand(mE,vE,m,n).*B;
b=uniform_rand(mb,vb,m,n).*B;

u=uniform_rand(mU,vU,m,1);
Beta=uniform_rand(mB,vB,m,1);
G=uniform_rand(mG,vG,n,1).*vectG';
g=uniform_rand(mg,vg,m,1);
mu_a=uniform_rand(mmA,vmA,n,1);
mu_p=uniform_rand(mmP,vmP,m,1);
w=uniform_rand(mw,vw,m,1);
phi=uniform_rand(mphi,vphi,m,1);

%Do not change the following line. Used to create globally accessible
%parameters
network_metadata = create_metadata(B, e, mu_p, mu_a, c, k, b, u, w, Beta, G, g, phi,bloomSpan,breakValue,flightSpan,breakValueF) ;
disp(G); disp(g) 

%Give initial variable states from uniform distribution 
mz=1.0; vz=1e-1; 
yzero=uniform_rand(mz,vz,2*m+n,1);

initial_plants=yzero(1:m);
initial_nectar=yzero(m+1:2*m);
initial_animals=yzero(2*m+1:2*m+n);
initial_alphas=B;

%Normalization an packing. You should not need to touch these two lines
initial_alphas=initial_alphas*diag(sum(initial_alphas).^(-1));
initial_alphas=initial_alphas(network_metadata.nz_pos) ;

%Integrate the model
initial_state=full([initial_plants;initial_nectar;initial_animals;initial_alphas]);
tspan = linspace(0,6000,6001);
options = odeset('JPattern', J_pattern) ;
[t y]=ode15s(@polinizacion_rhs_phenology,tspan,initial_state, options) ;

yi = y(1,:)';
yf = y(end,:)';

%Retrive variables' values across time
[pI, nectari, ai, alphasi] = unpack(yi, network_metadata);
[pf, nectarf, af, alphasf] = unpack(yf, network_metadata);
[plants, nectar, animals, avgAlphasf] = unpack2(y,m,n,network_metadata);
T=t;

%plot simulation output as time series
figure
subplot (3,1,1)
plot(T,plants)
subplot (3,1,2)
plot(T,nectar)
subplot (3,1,3)
plot(T,animals)

%Eliminate potential error 
assertP(all(all(alphasf>-1e-5 & alphasf<1.00001)));
assertP( all( abs(sum(alphasf) - 1.0)<1e-3 ) ) ;

assertP(all(all(avgAlphasf>-1e-5 & avgAlphasf<1.00001)));
%assertP( all( abs(sum(avgAlphasf) - 1.0)<1e-5 ) ) ;


