function [logfinal,final,obfngeo,obfnslip,objfnlikeli] = posteriorTohoku(model,others)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [logfinal,final] = posteriorTohoku(model,data,coord,disct,disp1,plot1)
% this function calculates the posterior probability of the model given the
% data. 
%
% outputs are: 
%   final = posterior probalility of the model given the data
%   logfinal = logarithm of the above posterior probability
%
% inputs are: 
%   model = the present state of model 
%           The first seven parameters are the geometrical parameters
%           The eighth and ninth parameters are hyperparameters
%           The tenth and the rest are the slip parameters    
%   data = geodetic data used to constrain the model 
%   los = line of sight of the geodetic data 
%   coord = coordinates of the geodetic datapoints 
%   LatEQ, LonEQ, depthEQ = previous earthquake data used to define the 
%                            prior probablity of the model 
%   disct = number of discretizations along fault dip
%
% In this function, we use the Least-square estimation to estimate the slip
% at the fault patches, and using the estimated slip, the posterior is cal-
% culated.
%
% written by : Rishabh Dutta, 22 Oct, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disct = 8; 
disp1 = others.disp1; 
plot1 = others.plot1;
surface = others.surface;
LatEQ = others.LatEQ;   % for geo prior
LonEQ = others.LonEQ;   % for geo prior
depthEQ = others.depthEQ;    % for geo prior
%greens = others.greens; 
subdisp = others.subdisp; 
subloc = others.subloc; 
weight = 0;     
LB = others.LB; 
UB = others.UB; 

model = reshape(model,1,length(model));
t = find(model <=UB & model >= LB); 

final = 0; 
logfinal = -inf; 

if length(t) ~= length(model)
    %fprintf('not within bounds\n')
    return 
end

% set number of discretization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(disct) == 1
    disct = 8; 
end
 
model1 = model(1:7); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% make the triangular element for the fault %%%%%%%%%%%%%%%%%%%%%%%%
 b1 = cputime;
 [p,q,r,trired,ang] =makemesh_full_inv(surface,model1,disct);
 % [p,q,r,trired,ang] =makemesh_full_inv(surface,bestgeo,disct);
 b2 = cputime;
 if ang == 0 
     return
 end
 
 if disp1 == 1
     fprintf('time taken to create the triangular mesh = %d\n',b2-b1);
 end
if plot1 == 1
    figure(1); trimesh(trired,p,q,r); axis equal; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% set the slip values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numpars = size(trired,1); 
slipall = reshape(model(11:end),2*numpars,1); 

dipslip =  slipall(1:numpars); 
strikeslip = slipall(numpars+1:end); 

if plot1 == 1
    figure(2); 
    subplot(121);patch(p(trired)',q(trired)',r(trired)',dipslip'); axis equal;
    subplot(122);patch(p(trired)',q(trired)',r(trired)',strikeslip'); axis equal;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define the hyperparameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
musq = model(8);            % for geometrical prior
rhosq = model(9);           % for slip prior
sigmasq = model(10);        % for likelihood 

alphasq = sigmasq/musq; 
betasq = sigmasq/rhosq; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the geometrical prior probabilty %%%%%%%%%%%%%%%%%%%%%%%%%
 b3 = cputime;
 model2 = model(1:8);
 [logprior1,prior1,obfngeo] = priorTohoku(model2,p,q,r,disct,LatEQ,LonEQ,depthEQ); 
 b4 = cputime; 
 if disp1 == 1
     fprintf('time taken for geometrical prior computation = %d\n',b4-b3);
     fprintf('Geometrical prior obj func: %d\n',obfngeo);
 end 
 if obfngeo == inf
     return
 end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the slip prior probability %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model3 = [reshape(model(1:7),1,7) model(9)];
b5 = cputime; 
[logprior2,prior2,obfnslip] = slipprior(model3,trired,p,q,r,dipslip,strikeslip);
b6 = cputime; 
if disp1 == 1
    fprintf('time taken for slip prior computation = %d\n',b6-b5);
    fprintf('Slip prior obj func: %d\n',obfnslip);
end 
if obfnslip == inf 
    return 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate greens function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b7 = cputime;
[greens,obsdata] = grn_func(subloc,subdisp,trired,p,q,r);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calculate the likelihood function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model4 = [dipslip; strikeslip]; 
if isempty(weight) == 1
    weight = 0; 
end

GPSdata = subdisp; 
GPScoord = subloc; 
[loglikelihood,likelihood,objfnlikeli] = likelinew(model4,greens,obsdata,...
                                            sigmasq,GPScoord,GPSdata,plot1); 
b8 = cputime; 
if disp1 == 1
    fprintf('time taken to compute likelihood = %d\n',b8-b7);
    fprintf('data obj func: %d\n',objfnlikeli);
end 
if objfnlikeli == inf
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

final = exp(-1/sigmasq*(objfnlikeli  + alphasq* obfngeo+ betasq *obfnslip));
logfinal = -1/sigmasq*(objfnlikeli  +alphasq* obfngeo + betasq * obfnslip);
b9 = cputime; 
if disp1 == 1
    fprintf('posterior value = %d\n log of posterior = %d\n',final,logfinal);
    fprintf('total time taken = %d\n',b9-b1); 
    fprintf('objfnlikeli = %f; \nbetasq*obfnslip= %f\n alpha*obfngeo = %f\n',...
        objfnlikeli,betasq * obfnslip,alphasq* obfngeo);
end 

end
    


















