function output = prior_samples(opt)

% initialisation : stage 0 
rng('shuffle');
N = opt.N;                         % Number of markov chains
Neff = opt.Neff;                   % chain length 
target = opt.target;               % target distribution
LB = opt.LB;                       % lower bound
UB = opt.UB;                       % upper bound
subdisp = opt.subdisp;
subloc = opt.subloc;

dims = length(LB);                 % dimensions
if size(LB,2)~=dims                % bounds should be put as row vectors
    fprintf('Error: Put bounds as row vector\n'); 
end
if sign(target(LB)) == 1           % target should be log of posterior
    fprintf('target disctribution has to be log posterior\n')
    return
end

% at stage 0; N number of samples are generated from the prior PDF 
diffbnd = UB-LB;
diffbndN = repmat(diffbnd,N,1);
LBN = repmat(LB,N,1);
smpstag0 = LBN + rand(size(diffbndN)).* diffbndN;   % input samples from PDF 

fprintf('loading the prior\n')

LBp = LB(1:7); UBp = UB(1:7);
diffbnd = UBp-LBp;
diffbndN = repmat(diffbnd,N,1);
LBN = repmat(LBp,N,1);
geoparsall = LBN + rand(size(diffbndN)).* diffbndN;
numpar = dims-10; 
slipsall = zeros(N,numpar);

parfor i = 1: N 
    bestgeo = geoparsall(i,:); 
    [slipvec] = lin_inv(bestgeo,subdisp,subloc); 
    slipsall(i,:) = slipvec; 
end
fprintf('prior run finished\n')
hyperall = [repmat(LB(8),N,1) repmat(LB(9),N,1) repmat(LB(10),N,1)]; 
smplstage0 = [geoparsall hyperall slipsall];  

smplstage0 = smplstage0(1:N,:); 
smplstage0(find(smplstage0 > -1e-5 & smplstage0 < 1e-5)) = 0; 
smplstage0(:,8) = LB(8);
smplstage0(:,9) = LB(9);
smplstage0(:,10) = LB(10);
clear geoparsall slipsall hyperall

save('samples_taken1_1','smplstage0')

logpst = zeros(N,1);
parfor i = 1:N          % calculate posterior of samples
    smp0 = smplstage0(i,:); 
    logpost = target(smp0); 
    logpst(i,1) = logpost; 
end

% find beta for stage 0 
post = exp(logpst);                
beta = rand*.1;                    % random beta
wght = post.^beta;
covwght = std(wght)/mean(wght);    % coefficient of variation of posteriors
refcov = covwght; 
fracrefcov = .01* refcov; 
run = 1; 
while covwght> (refcov+fracrefcov) || covwght< (refcov-fracrefcov) % coeff of var should be close to 1
    beta = rand; 
    wght = post.^beta; 
    covwght = std(wght)/mean(wght); 
    run =run +1; 
    if mod(run,5000) ==1 
        fprintf('running more than 5000 times; run again \n');
    end        
end

probwght = wght/sum(wght);         % probabilty of the samples
normwght = wght/max(wght);         % norm weight of samples

% resampling at the stage 0 
inind = 1:N;
qwht = probwght;
outind = deterministicR(inind,qwht);

%     i = 0; indwght = [];
%     while i < N
%         ind = ceil(rand(N-i,1)*N);
%         resampwght = probwght(ind');
%         comparewght = resampwght - rand(N-i,1);
%         ind2 = find(comparewght>0);
%         indwght = [indwght ind(ind2)'];
%         i = i+length(ind2);
%     end
% resampled samples
resmpl = smplstage0(outind,:);

save('samples_taken1_2','smplstage0','resmpl','outind','logpst','beta')
output.resmpl = resmpl; 
