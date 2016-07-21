% run for the linear inversions 

clear
%addpath additional_scripts
load parameters
disct =8;
surface = surf_trace; 
load EQ_data1
load GPS_subsampled

bestgeo = [0.0571    7.1445   58.0000    0.8426   0.0912    1.2809    -0.1008]; 
%bestgeo = [0.0571    7.1445   58.0000    1  .5    1.2809    0.1008]; 
[p,q,r,trired,ang] =makemesh_full_inv(surface,bestgeo,disct);
%trimesh(trired,p,q,r); axis equal

[grn1,obsdata] = grn_func(subloc,subdisp,trired,p,q,r);


lowslip1 = zeros(1,size(trired,1));
lowslip2 = -5*ones(1,size(trired,1));
maxslip1 = 50*ones(1,size(trired,1));
maxslip2 = 3*ones(1,size(trired,1));
LBslip = [lowslip1 lowslip2]; UBslip = [maxslip1 maxslip2];
LB = [.05 7 58 0 -0.3 0 -0.3 3.5e4*2 1/2e7*2 2  LBslip]'; 
UB = [.07 10 58 5 0.3 5 0.3 3.5e4*2 1/2e7*2 2 UBslip]';
LB(5)

edges= [];
for i = 1:10
    if i == 1
        edges = [1:disct*2];
    elseif i<10
        edges = [edges i*disct*2-1 i*disct*2];
        %edges = [edges (i-1)*disct*2+1 (i-1)*disct*2+2 i*disct*2-1 i*disct*2];
    elseif i==10
        edges = [edges (i-1)*disct*2+1:i*disct*2];
    end
end
    
LB(10+size(trired,1)+edges) = 0; 
UB(10+size(trired,1)+edges) = 0;

% Define posterior as a function of model ................
others.disct = 8; 
others.disp1 = 0; 
others.plot1 = 0; 
others.weight = 0; 
others.LatEQ = LatEQ2; 
others.LonEQ = LonEQ2; 
others.depthEQ = depth2; 
others.data = []; 
others.coord = [];
others.LB = LB'; 
others.UB = UB'; 
others.obsdata = obsdata; 
others.subdisp = subdisp;
others.subloc = subloc; 
others.surface = surface; 
%others.greens = grn1; 
%others.trired = trired; 
%others.p = p; 
%others.q = q; 
%others.r = r; 
posterior = @(x) posteriorTohoku(x,others);

% Define posterior as a function of model for plot and time
other1.disct = 8; 
other1.disp1 = 1; 
other1.plot1 = 1; 
other1.weight = 0; 
other1.LatEQ = LatEQ2; 
other1.LonEQ = LonEQ2; 
other1.depthEQ = depth2; 
other1.data = []; 
other1.coord = [];
other1.LB = LB'; 
other1.UB = UB';
other1.obsdata = obsdata;
other1.subdisp = subdisp;
other1.subloc = subloc; 
other1.surface = surface; 
%other1.greens = grn1; 
%other1.trired = trired; 
%other1.p = p; 
%other1.q = q; 
%other1.r = r;
plotpost = @(x) posteriorTohoku(x,other1);

%%
opt.N  = 5*2; 
opt.Neff = 6
opt.target = posterior; 
opt.LB = LB'; 
opt.UB = UB'; 
opt.subdisp = subdisp; 
opt.subloc = subloc; 
fprintf('done2\n')
 
tic; outsmpl = prior_samples(opt);  b = toc; 



