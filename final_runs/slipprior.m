function [logfinal,final, prior1] = slipprior(model,trired,p,q,r,dipslip,strikeslip)
% function [logfinal,final] = slipprior(model,trired,p,q,r,dipslip,strikeslip)
% this function calculates the prior propability of the present state of
% slip model using the laplacian.
% outputs are : 
%   final= the prior value 
%   logfinal = log of prior
% inputs are : 
%   model = present state of model
%   dipslip, strikeslip = dipslip and strikeslip values
% written by : Rishabh Dutta
% Edited: 26 Oct, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

slipall = [dipslip(:); strikeslip(:)];
betasq = model(8);

% get the laplacian %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
laplac1 = laplacian(trired,p,q,r); 
laplac = [laplac1 zeros(size(laplac1)); zeros(size(laplac1)) 1*laplac1]; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lapslip = laplac*slipall;
final2 = (-1/(2*betasq))*(lapslip'*lapslip);
wgtslip = norm(slipall)^-2;
final3 = wgtslip * final2; 

final = betasq^(-.5)*exp(final3); 
logfinal = -.5*log(betasq)+ final3; 

prior1 = (lapslip'*lapslip)*wgtslip;
end


