function output = ATMIP_remaining2(opt,beta,samplestage,stage,num,covsmpl) 

fprintf('inside the script ATMIP_remaining 2\n')
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


% Adaptive metropolis algorithm from N chains
mhsmpl = zeros(N,dims);
mhpost = zeros(1,N);
acc_rate = zeros(1,N);

fprintf('At stage %d; beta = %f; Started metropolis chains\n',stage,beta);

parfor i =1:N
    start = samplestage(i,:);
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
save([dirname '/' subdir '/sol_3rd.mat'],'opt','beta','post','samplestage','stage','covsmpl');

fileID = fopen([dirname '/' subdir '/readme.txt'],'a');
fprintf(fileID,'samples taken with next Neff = %d \n stage = %d, beta = %f\n',Neff,stage,beta); 
fclose(fileID);

% final results 
output.samples = mhsmpl; 
output.posterior = mhpost; 
output.stages = stage; 