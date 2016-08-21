function [logfinal, final, prior1] = priorTohoku(model,p,q,r,disct,LatEQ,LonEQ,depthEQ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [logfinal, final] = priorTohoku(model,p,q,r,disct,LatEQ,LonEQ,depthEQ)
% this function calculates the prior propability of the present state of
% model using the previous earthquakes data.
% outputs are : 
%   final= the prior value 
%   logfinal = log of prior
% inputs are : 
%   model = present state of model
%   disct = number of discretization of fault along dip
%   LatEQ, LonEQ, depthEQ = all previous earthquake data
% written by : Rishabh Dutta
% Edited: 26 Oct, 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load parameters
% disct =15; 
% coefd = model(1:3);
% strike = model(4);
% maxdep = model(5);

alphasq = model(8);
% [p,q,r,trired] = makemesh2(surf_trace,coefd,strike,maxdep,disct);
nofault = length(p); 
disct1 = disct+1; 

% addpath /home/duttar/Desktop/insar/japan/Tohoku/Tohoku/util/
% addpath /home/duttar/Desktop/insar/japan/Tohoku/Tohoku/data/
% load EQdata
% load EQ_ind
% 
% LongitudeE = LongitudeE(ind);
% LatitudeN = LatitudeN(ind);
% JMADepthkm = JMADepthkm(ind);
% 
% ref_coord = [140 38]; 
% [xy_EQ]=llh2localxy([LongitudeE';LatitudeN'],ref_coord);
% LatEQ = -xy_EQ(:,1); 
% LonEQ = xy_EQ(:,2) ; 
% depthEQ = -JMADepthkm ;

noEQ = length(LatEQ);
distanc = zeros(noEQ,1);

parfor i=1:noEQ
    xpoint = LonEQ(i); xpoint = xpoint + randn; 
    ypoint = LatEQ(i); ypoint = ypoint + randn; 
    zpoint = depthEQ(i); zpoint = zpoint + 3*randn; 
    % find the closest points on the fault      
    distall = sqrt((r-zpoint*ones(nofault,1)).^2+...
        (q-ypoint*ones(nofault,1)).^2+...
        (p-xpoint*ones(nofault,1)).^2);   
    sortdis = sort(distall);     
    ind = zeros(length(nofault),1);
    for j = 1:nofault
        ind(j) = find(distall == sortdis(j)); 
    end   
    ind2 = ind;
    colall = ceil(ind./disct1);     
    if colall(1) == colall(2)
        n=3;
        while colall(3) == colall(2)
             colall(3) = colall(n+1);
            ind2(3) = ind2(n+1);
            n=n+1; 
        end
    else
        n=3; 
        while (colall(3) ~= colall(1)) &&  (colall(3) ~= colall(2))            
                colall(3) = colall(n+1);
                ind2(3) = ind2(n+1);
                n=n+1;       
        end
    end
    
    % so the three points of the plane are 
    planept1 = [p(ind2(1)) q(ind2(1)) r(ind2(1))];
    planept2 = [p(ind2(2)) q(ind2(2)) r(ind2(2))];
    planept3 = [p(ind2(3)) q(ind2(3)) r(ind2(3))];
    
    normal = cross(planept1-planept2,planept1-planept3);
    A = normal(1); B = normal(2); C = normal(3); 
    D = -dot(normal,planept2);
    distanc(i) = abs(A*xpoint+B*ypoint+C*zpoint+D)/sqrt(A^2+B^2+C^2); 
    
%    dist1 = sqrt((xpoint-planept1(1))^2+(ypoint-planept1(2))^2+(zpoint-planept1(3))^2);
%    dist2 = sqrt((xpoint-planept2(1))^2+(ypoint-planept2(2))^2+(zpoint-planept2(3))^2);
%    dist3 = sqrt((xpoint-planept3(1))^2+(ypoint-planept3(2))^2+(zpoint-planept3(3))^2);
    
%     normal = cross(planept1-planept2, planept1-planept3);
%     syms xy yz zx;
%     planept = [xy yz zx];
%     %planefunction = dot(normal, planept-planept1);
%     realdot = @(u, v) u*transpose(v);
%     planefunction = realdot(planept-planept1,normal);
%     planefunc = matlabFunction(planefunction);
%     
%     %coefficients in the plane equation 
%     D = planefunc(0,0,0);
%     A = planefunc(1,0,0) - D; 
%     B = planefunc(0,1,0) - D;
%     C = planefunc(0,0,1) - D;
%     smallestdist = (A*xpoint+B*ypoint+C*zpoint+D)/sqrt(A^2+B^2+C^2);    
%     distanc(i) = (dist1+dist2+dist3)/3; 
end            
    
final2 = sqrt(sum(distanc.^2));     
    
final = alphasq^(-.5)*exp(-1/(2*alphasq)* final2.^2) ; 
logfinal = -.5* log(alphasq) - 1/(2*alphasq) * final2.^2;
    
prior1 = final2.^2; 
% defined from looking at the geometry of the fault

% lowlimit =  [.02 7 55 -10 -1.5 -10 -1.5 0];
% uplimit = [.09 15 70 10 1.5 10 1.5  1e12];
% 
% t=find(model<uplimit & model>lowlimit); 
% if length(t) ~= length(model)
%     final = 0; 
%     logfinal = -inf;
% end 
end
    
    
    
    






