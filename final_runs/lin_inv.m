function [slipvec,trired,p,q,r] = lin_inv(model,subdisp,subloc)
%  slipvec = lin_inv(trired,p,q,r,subdisp,subloc)

load parameters
disct =8; 
surface = surf_trace; 
bestgeo  = model(1:7); 

[p,q,r,trired,ang] =makemesh_full_inv(surface,bestgeo,disct);

lowslip1 = zeros(1,size(trired,1));
lowslip2 = -5*ones(1,size(trired,1));
maxslip1 = 50*ones(1,size(trired,1));
maxslip2 = 3*ones(1,size(trired,1));
LBslip = [lowslip1 lowslip2]; UBslip = [maxslip1 maxslip2];
LB = [LBslip]'; 
UB = [UBslip]';

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
LB(size(trired,1)+edges) = 0; 
UB(size(trired,1)+edges) = 0;


GPSdata = subdisp; 
GPScoord = subloc; 

numpars = 2*size(trired,1);
numdata = numel(GPSdata); 

tic
[greens,datavector] = grn_func(GPScoord,GPSdata,trired,p,q,r);
ltime = toc;
greens1 = greens;

musq = 30;
laplac1 = laplacian(trired,p,q,r); 
laplac = musq*[laplac1 zeros(size(laplac1)); zeros(size(laplac1)) laplac1]; 

Amat = [greens1; laplac];  % (numdata+ 0.5*numpars)*numpars
Bmat = [datavector; zeros(numpars,1)];

tikhA = Amat'*Amat + .005^2*eye(size(Amat'*Amat));
options = optimoptions('lsqlin','display','off','Algorithm','active-set');
[slipvec,resnorm,residual] = lsqlin(tikhA,Amat'*Bmat,[],[],[],[],LB,UB,[],options);

