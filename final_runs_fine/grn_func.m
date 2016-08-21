function [greens,datavector] = grn_func(GPScoord,GPSdata,trired,p,q,r)
% [greens,datavector] = grn_func(GPScoord,GPSdata,trired,p,q,r)
% forms the greens function for the GPS data coordinates with the geometry
% defined by trired indices 
% outputs are greens function (Note: This matrix has to be combined with
% laplacian) and datavector
% inputs are the GPS data and coordinates and trired indices and fault
% nodes. 
% [greens,datavector] = grn_func(GPScoord,GPSdata,trired,p,q,r)
% written by Rishabh Dutta, 28 Sep 2015

% Add meade's library
%addpath /home/duttar/Desktop/softwares/triangular_dislocation/

numdata = numel(GPSdata); 
numpars = 2*size(trired,1);    % we have total of 2*numpars slip parameters

p_tri = p(trired); q_tri = q(trired); r_tri = -r(trired); 
xcoord = GPScoord(:,1); 
ycoord = GPScoord(:,2);
zcoord = -GPScoord(:,3); 

greens1= zeros(numdata,numpars/2); % this refers to green's function due to 
                                 % thrust-slip component
greens2= zeros(numdata,numpars/2); % this refers to green's function due to 
                                 % strike-slip component                                
parfor i=1: numpars/2
    % element coordinates:
    xparco = p_tri(i,:); yparco = q_tri(i,:); zparco = r_tri(i,:); 
    Uall1 = CalcTriDisps(xcoord,ycoord,zcoord,xparco,...
        yparco,zparco,.28,0,0,-1);
    dataunit1 = [Uall1.x(:) Uall1.y(:) -Uall1.z(:)]; % set as same config as GPSdata
    greens1(:,i) = dataunit1(:);
    
    Uall2 = CalcTriDisps(xcoord,ycoord,-zcoord,xparco,...
        yparco,-zparco,.28,1/sqrt(2),0,1/sqrt(2));
    dataunit2 = [Uall2.x(:) Uall2.y(:) -Uall2.z(:)]; % set as same config as GPSdata
    greens2(:,i) = dataunit2(:); 
end

greens = [greens1 greens2]; 
datavector = reshape(GPSdata(:),numdata,1); 


