% Adaptive Metropolis algorithm 
% scales the covariance matrix according to the acceptance rate 
% cov of proposal = (a+bR)*sigma ; R = acceptance rate 
% returns the last sample of the chain
%
% syntax: [G,GP] = AMH(X,target,covariance,mrun)
%
% INPUTS: 
% - X (starting model)
% - target (target distribution)
% - covariance (covariance of the proposal distribution)
% - mrun (number of samples)
% - beta (samples posterior^beta )
%
% OUTPUTS: 
% - G (final sample of the chain)
% - GP (posterior value of the final sample)
%
% written by: Rishabh Dutta, 12 Mar 2016
% (Don't forget to acknowledge)

function [G,GP,avg_acc] = AMH(X,target,covariance,mrun,beta,LB,UB)

Dims = length(X);       % dimensions of the model
logpdf = target(X);     % log of posterior value of starting model
V = covariance;         % covariance of proposal distribution
best_X = X;             % starting model
best_P = logpdf*beta;        
P0 = logpdf*beta; 
eps = -1; 
pcom = 0; 
a = 1/180;        % for scaling proposal covariance
b = 179/180;        % for scaling proposal covariance

sameind = find(LB==UB); 

% Define optimal acceptance rate
switch Dims
    case 1
        optimal=0.441;
    case 2
        optimal=0.352;
    case 3
        optimal=0.316;
    case 4
        optimal=0.285;
    case 5
        optimal=0.275;
    case 6
        optimal=0.273;
    case 7
        optimal=0.270;
    case 8
        optimal=0.268;
    case 9
        optimal=0.267;
    case 10
        optimal=0.266;
    case 11
        optimal=0.265;
    case 12
        optimal=0.264;
    otherwise
        optimal=0.255;
end

% Set initial scaling factor
s = a+ b*optimal;

% metropolis chain
mark = cputime; mark_0 = mark;

U = log(rand(1,mrun));    
TH = zeros(Dims,mrun);
THP = zeros(1,mrun);
avg_acc = 0;
factor = zeros(1,mrun);
for i = 1:mrun
    X_new = mvnrnd(reshape(X,1,length(X)),s^2*V);  % new sample
    X_new(sameind) = LB(sameind); 
    ind1 = find(X_new < LB); 
    diff1 = LB(ind1) - X_new(ind1); 
    X_new(ind1) = LB(ind1) + diff1; 
    
    ind2 = find(X_new > UB); 
    diff2 = X_new(ind2) - UB(ind2); 
    X_new(ind2) = UB(ind2) - diff2; 
    
    P_new = beta*target(X_new);                         % posterior value at new sample
    if P_new > best_P                              % accept sample if Pnew>P0
        best_X=X_new;
        best_P = P_new;
        X = X_new;
        P0 = P_new;
        acc_rate=1;                                % acceptance rate
    else
        rho = P_new - P0;
        acc_rate=exp(min(rho,0));
        if U(i)<= rho                              % accept X_new with a probability
            X = X_new;                             
            P0 = P_new;
        end
    end
    
    TH(:,i) = X;
    THP(i) = P0;
    factor(i) = s^2;            
    avg_acc = avg_acc*(i-1)/i + acc_rate/i;        % calculate acceptance rate
    s = a + b*avg_acc;                             % scale for proposal covariance updated
end

% results
G = TH(:,end)';         % final sample 
GP = THP(end)/beta;          % posterior value


















