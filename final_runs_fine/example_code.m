model1 = [ 0.0562    9  60 .1 -.05 -1.5 -.05 1e-5 5 1e2];
model2 = 20 * rand(1,336);
model3 = 5 * randn(1,336);
model = [model1 model2 model3];

addpath ../extras/
disct =8;
[slip1,slip2] = make_slip(trired,p,q,r,disct+1);
model = [model1 slip1(:,1)' slip2(:,1)'];


[p,q,r,trired,ang] =makemesh_full_inv(surface,model1,disct);
figure(3); trimesh(trired,p,q,r); axis equal; 
tic ;[logval,val]= posteriorTohoku(model,[],[],8,1,1,LatEQ,LonEQ,depthEQ,0); toc

%% check the displacement 
clear; clc
disct  = 8; 
model1 = [ 0.0793    8.87  58 .1 -.05 -1.5 -.05 1e-5 5 1e2];
load /home/duttar/Desktop/insar/japan/Tohoku/Tohoku/modeling/Nonlinnew/new_pars/full_inversion/full_bayesian/data/GPSall.mat
load /home/duttar/Desktop/insar/japan/Tohoku/Tohoku/modeling/Nonlinnew/new_pars/full_inversion/full_bayesian/data/parameters.mat
load /home/duttar/Desktop/insar/japan/Tohoku/Tohoku/modeling/Nonlinnew/new_pars/full_inversion/full_bayesian/data/varGPS.mat
surface = surf_trace; 
[p,q,r,trired,ang] =makemesh_full_inv(surface,model1,disct);
[slip1,slip2] = make_slip(trired,p,q,r,disct+1);
model = [model1 slip1(:,1)' slip2(:,1)'];
dipslip = slip1(:,1); strikeslip = slip2(:,2); 
addpath /home/duttar/Desktop/softwares/triangular_dislocation/

p_tri = p(trired); q_tri = q(trired); r_tri = -r(trired); 
xcoord = GPScoord(:,1); 
ycoord = GPScoord(:,2);
zcoord = -GPScoord(:,3); 

parfor i =1: size(trired,1)
    xparco = p_tri(i,:); yparco = q_tri(i,:); zparco = r_tri(i,:); 
    sslip = strikeslip(i); dslip = dipslip(i);
    Uall = CalcTriDisps(xcoord,ycoord,zcoord,xparco,...
        yparco,zparco,.26,sslip,0,-dslip);
    dataunit = [Uall.x(:) Uall.y(:) -Uall.z(:)];
    preddata1(:,i) = dataunit(:); 
end

preddata2 = preddata1';
preddata = sum(preddata2);

xpred = preddata(1:1234); ypred = preddata(1235:2468); 
zpred = preddata(2469:3702); 

subplot(211); quiver(GPScoord(:,1),GPScoord(:,2),xpred',ypred')
axis equal; axis([-1230 1400 -700 600])
subplot(212); scatter(GPScoord(:,1),GPScoord(:,2),10,zpred')
axis equal; axis([-1230 1400 -700 600])



%%

configCluster

smc = parcluster('smc'); 
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('st_0_0')
ClusterInfo.setQueueName('short')
ClusterInfo.setWallTime('240')

job = batch(smc, 'change_stagerun', 'pool',449, 'AttachedFiles', 'additional_scripts/');
job = batch(smc, 'nexthalf_run', 'pool',449, 'AttachedFiles', 'additional_scripts/');

job = batch(smc, 'stage0_priors', 'pool',249, 'AttachedFiles', 'additional_scripts/');


job = batch(smc,'run_tmcmc2', 'pool',99, 'AttachedFiles', 'additional_scripts/');


job = batch(smc,'firsthalf_run', 'pool',249, 'AttachedFiles', 'additional_scripts/');


job = batch(smc, 'tmcmc_runagain834', 'pool',499, 'AttachedFiles', 'additional_scripts/');



noor1 = parcluster('noor1'); 
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('run_nogeo')
ClusterInfo.setQueueName('rh6_q2hr')
ClusterInfo.setWallTime('2:00')

job = batch(noor1, 'run_tmcmc2', 'pool',255, 'AttachedFiles', 'additional_scripts/');


noor1 = parcluster('noor1'); 
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('run_nogeo')
ClusterInfo.setQueueName('rh6_q2hr')
ClusterInfo.setWallTime('2:00')

job = batch(noor1, 'run_tmcmc2', 'pool',255, 'AttachedFiles', 'additional_scripts/');

noor2 = parcluster('noor2'); 
ClusterInfo.setEmailAddress('rishabh.dutta@kaust.edu.sa')
ClusterInfo.setJobName('run_par')
ClusterInfo.setQueueName('defaultq')
ClusterInfo.setWallTime('240')

job = batch(noor2, 'run_tmcmc2', 'pool',511, 'AttachedFiles', 'additional_scripts/');


