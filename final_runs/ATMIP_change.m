function output = ATMIP_change(opt,beta,post,samplestage,stage,num) 

fprintf('inside the script ATMIP_remaining \n')
% initialisation : stage 0 
rng('shuffle');
N = opt.N;                         % Number of markov chains
Neff = opt.Neff;
Neff = 40;                     % chain length 
target = opt.target;               % target distribution
LB = opt.LB;                       % lower bound
UB = opt.UB;                       % upper bound
dims = length(LB);                 % dimensions
if size(LB,2)~=dims                % bounds should be put as row vectors
    fprintf('Error: Put bounds as row vector\n'); 
end
if sign(target(LB)) == 1           % target should be log of posterior
    fprintf('target disctribution has to be log posterior\n')
    return
end

% save the data for stage 0 
randno = num; 
randno = num2str(randno);
dirname = 'runagain'; 
dirname = strcat(dirname,randno); 
[s,mess,messid] = mkdir(dirname);

fprintf(mess)
fprintf('\n')
 
% beta value goes from 0 to 1

beta1 = beta; 
if stage ~= 0 
    beta = beta1 + .1;
    while beta < beta1   % random beta
        beta = rand; 
    end
    wght = post.^(beta-beta1);
    covwght = std(wght)/mean(wght);    % coefficient of variation of posteriors
    refcov = 1; 
    fracrefcov = .01* refcov;
    n=1; 
    while covwght> (refcov+fracrefcov) || covwght< (refcov-fracrefcov) % coeff of var should be close to 1
        beta = rand;
        while beta < beta1   % random beta
            beta = rand;
        end
        wght = post.^(beta-beta1); 
        covwght = std(wght)/mean(wght); 
        if n > 5000 && covwght<refcov 
            beta = 1; 
            break
        end
        n=n+1; 
    end
    
    refcov
    covwght
    beta = min(1,beta);   % beta at stage m+1
    beta1 = beta;         % beta at stage m 
    
    probwght = wght/sum(wght);         % probabilty of the samples
    normwght = wght/max(wght);         % norm weight of samples
    
    % resampling
    
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
    resmpl = samplestage(outind,:);
    
    % define the covariance for the proposal distribution used in AMH sampling
    meansmpl = sum(repmat(probwght,1,dims).* samplestage);
    
    covsmpl = zeros(dims,dims);
    tic; 
    for i = 1:N
        par = samplestage(i,:);
        smpldiff = par - meansmpl;
        covint = probwght(i)* (smpldiff'*smpldiff);
        covsmpl = covsmpl+covint;
        b= toc; 
        if b >100 && b<110
            fprintf('creating cov for MH, taking time \n')
        end
    end
    
    % covariance of the samples
   
    
    % Adaptive metropolis algorithm from N chains
    mhsmpl = zeros(N,dims);
    mhpost = zeros(1,N);
    acc_rate = zeros(1,N);

    fprintf('At stage %d; beta = %f; Started metropolis chains\n',stage,beta);
    
    parfor i =1:N        
        start = resmpl(i,:);
        [G,GP,acc_r] = AMH(start,target,covsmpl,Neff,beta,LB,UB);
        mhsmpl(i,:) = G;
        mhpost(i) = GP;
        acc_rate(i) = acc_r;    
    end
    nomoves = length(find(acc_rate==0));
    fprintf('MH didnt change for %d chains at stage %d\n',nomoves,stage);
    
    post = exp(mhpost'); 
    samplestage = mhsmpl; 
    beta1 = beta; 
    
    % save the intermediate data
    subdir = strcat('stage',num2str(stage)); 
    mkdir(dirname,subdir); 
    save([dirname '/' subdir '/sol_1st.mat'],'opt','beta','post','samplestage','covwght','stage','covsmpl','num');
    
        
end

fileID = fopen([dirname '/' subdir '/readme.txt'],'w');
fprintf(fileID,'samples taken with first Neff = %d \n stage = %d, beta = %f\n',Neff,stage,beta); 
fclose(fileID);

% final results 
output.samples = mhsmpl; 
output.posterior = mhpost; 
output.stages = stage; 