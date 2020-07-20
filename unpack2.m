function [plants, nectar, animals, avgAlphasAssigned] = unpack2(y,m,n, network_metadata)

m = network_metadata.plant_qty ;
n = network_metadata.animal_qty ;

plants=y(:,1:m) ;
nectar=y(:,m+1:2*m);
animals=y(:,2*m+1: 2*m+n);

avgAlphas=mean(y(2000:3000,:),1)';%I think the alphas change quickly...therefore 2*n+m+1:end
%when we look at the end of the sims, we're seeing more specialized than we
%should be. SO I'm looking at the avg alphas across the last 1/3 of the sim

avgAlphasAssigned=sparse( m, n ) ;
avgAlphasAssigned(network_metadata.nz_pos) = avgAlphas(2*m+n+1:end,1) ;

end
