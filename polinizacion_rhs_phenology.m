%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code created for integrating phenological modifiers into plant and
% pollinator behavior and interactions in plant-pollinator networks
% using Valdovinos et al. 2013, 2016 model as basis
% Original Model Authors: Fernanda S. Valdovinos & Pablo Moisset de Espanes
% Modified to include phenology dynamics by Paul Glaum 
% Last Modification: July, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model can be switched to running with no phenology by switching commented
% lines (see below). 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dx = polinizacion_rhs_phenology(t,x)

global network_metadata indRemP indRemA

plant_qty  = network_metadata.plant_qty ;
animal_qty = network_metadata.animal_qty ;
nz_pos = network_metadata.nz_pos ;
e      = network_metadata.e ;
mu_p   = network_metadata.mu_p ;
mu_a   = network_metadata.mu_a ;
c      = network_metadata.c ;
b      = network_metadata.b ;
u      = network_metadata.u ;
w      = network_metadata.w ;
Beta   = network_metadata.Beta ;
G      = network_metadata.G ;
g      = network_metadata.g ;
phi    = network_metadata.phi ;
In     = network_metadata.In ;
bloomSpan=network_metadata.bloomSpan ;
breakValue=network_metadata.breakValue ;
flightSpan=network_metadata.flightSpan ;
breakValueF=network_metadata.breakValueF ;

[p N a Alpha] = unpack(x, network_metadata ) ;
p(indRemP)=0;% We force extinct species to zero to avoid resucitation due to stiff integrator
N(indRemP)=0;
a(indRemA)=0;

%% Phenology wave functions
m=plant_qty;%the # of plants
bloom=1/bloomSpan;%baseline wave freq, controls width of floral peak through bloomSpan
bloombreak=round(breakValue*2)/2;% # of peaks to skip b4 bloom starts again, using breakvalue
bfreq=(bloombreak/2)+.5;%how many periods to wait between blooms

d=1; %d=0 means complete overlap, d=1 means maximum separation possible between basal peaks
sepDist=((2*pi*(breakValue+1)/m))/d;%divide by the # of plants (m)for least complimentarity
separations=transpose((0:(m-1))*sepDist);
waveShiftSqAugment=bloombreak+1;

n=animal_qty; %# of animals 
flight=1/flightSpan;
flightbreak=round(breakValueF*2)/2;
ffreq=(flightbreak/2)+.5;

d=1; %d=0 means complete overlap, d=1 means maximum separation possible between basal peaks
sepDistF=((2*pi*(breakValueF+1)/n))/d;%divide by the # of plants (m)for least complimentarity
separationsF=transpose((0:(n-1))*sepDistF);
waveShiftSqAugmentF=flightbreak+1;

%Plant Phenology equation, Tpi...
Tpi=abs( (bloombreak+2)*sin((bloom*pi*t)-separations).*(1+square(((bloom/(2*bfreq))*pi*t)-(separations/waveShiftSqAugment),(25/bfreq)))/2 );
%Pollinator Phenology equation, Tai...
Tai=abs( (flightbreak+2)*sin((flight*pi*t)-separationsF).*(1+square(((flight/(2*ffreq))*pi*t)-(separationsF/waveShiftSqAugmentF),(25/ffreq)))/2 );

%% Model's specific computation begins here

sigma = diag(sparse(p)) * Alpha ; 
sigma = sigma * diag(sparse(1./(sum(sigma)+realmin))) ;

%no phenology
Gamma = g .* (1 - u'*p - w.*p + u.*p) ;% non sparse 

%pollination interactions, choose between active phenology and no phenology
%no phenology
%tmp = Alpha * diag(sparse( a_ph.*a ) ) ; %tmp nxm sparse
%With Phenology function
tmp = Alpha * diag(sparse( Tai.*a ) ) ; %tmp nxm sparse,
                        %  ^ here for bee phenology

dp = ( (Gamma .* sum(e .* sigma .* tmp, 2)) - mu_p) .* p;
tmp = (diag(sparse(N)) * tmp) .* b ;
da = sum(c .* tmp, 1)' - mu_a .* a ;  


%Rewards production, choose between active phenology and no phenology
%no phenology
%dN = Beta .* p - phi.*N - sum(tmp, 2) ; 
%With Phenology Function
dN = Tpi.* Beta .* p - phi.*N - sum(tmp, 2) ; 

% Adaptation starts here

%Fitness function
DH = diag(sparse(N)) * sparse(c.*b) ; %nxm sparse
DH(Alpha<0)=-DH(Alpha<0) ; %From Fernanda's GitHub

wavg = sum(Alpha.*DH) ; %Weights for average. nxm sparse

%This is the replicator equation, choose between active phenology and no phenology
dAlpha = Alpha.*DH - Alpha*diag(sparse(wavg)) ;
%no phenology
%dAlpha = dAlpha*diag(sparse( G )) ;      
%With Phenology Function
dAlpha = dAlpha*diag(sparse( Tai.*G )) ;    


%% Now pack the answer
dx = full([dp; dN; da; dAlpha(nz_pos)]) ;
