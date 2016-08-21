
function [laplac] = laplacian(trired,p,q,r)
% forms the laplacian for the tohoku fault geometry 
% Note: this creates laplacian for just one strikeslip or thrust-slip
% [laplac] = laplacian(trired,p,q,r)
% written by Rishabh Dutta, 28Sep 2015

% no. of patches 
npat = size(trired,1);
laplac = sparse(npat,npat); 

for i = 1: npat 
    % 3 corners of ith patch 
    indi = trired(i,:);
    centr_i = [mean(p(indi)) mean(q(indi)) mean(r(indi))];
    
    % now find the 3 triangles sharing the edges
    
    % 1st edge is in following patches: 
    firedge = mod(find(trired == indi(1)),npat); 
    firedge(firedge ==0) = npat; 
    
    % 2nd edge is in following patches: 
    secedge = mod(find(trired == indi(2)),npat);
    secedge(secedge ==0) = npat; 
    
    % 3rd edge is in following patches: 
    thiedge = mod(find(trired == indi(3)),npat); 
    thiedge(thiedge ==0) = npat; 
    
    % find the triangle sharing 1st and 2nd corners
    tri12 = []; 
    for j = 1: length(firedge)
        shr1 = secedge(find(secedge == firedge(j))); 
        if ~isempty(shr1) && shr1 ~= i
            tri12 = [tri12 shr1]; 
        end
    end
    
    % find the triangle sharing 2nd and 3rd corners
    tri23 = []; 
    for j = 1: length(secedge)
        shr1 = thiedge(find(thiedge == secedge(j))); 
        if ~isempty(shr1) && shr1 ~= i
            tri23 = [tri23 shr1]; 
        end
    end
    
    tri31 = []; 
    for j = 1: length(thiedge)
        shr1 = firedge(find(firedge == thiedge(j))); 
        if ~isempty(shr1) && shr1 ~= i
            tri31 = [tri31 shr1]; 
        end
    end
    
    tris = [tri12 tri23 tri31]; 
    
    if length(tris) == 3
        % center of 1st triangle 
        indvert = trired(tris(1),:);  
        centr_fir = [mean(p(indvert)) mean(q(indvert)) mean(r(indvert))];         
        distri1sq =  (centr_fir(1)-centr_i(1))^2 +(centr_fir(2)-centr_i(2))^2+...
            (centr_fir(3)-centr_i(3))^2; 
        % distri1sq = 3*(rms(centr_fir-centr_i).^2);
        
        % center of 2nd triangle
        indvert = trired(tris(2),:); 
        centr_sec = [mean(p(indvert)) mean(q(indvert)) mean(r(indvert))]; 
        distri2sq = 3*(rms(centr_sec-centr_i).^2);
        
        % center of 3rd triangle
        indvert = trired(tris(3),:); 
        centr_thi = [mean(p(indvert)) mean(q(indvert)) mean(r(indvert))]; 
        distri3sq = 3*(rms(centr_thi-centr_i).^2);
        
        laplac(i,tris(1)) = 1/distri1sq ; 
        laplac(i,tris(2)) = 1/distri2sq ;
        laplac(i,tris(3)) = 1/distri3sq ; 
        laplac(i,i) = -(1/distri1sq + 1/distri2sq + 1/distri3sq);
    
    elseif length(tris) == 2
        % center of 1st triangle 
        indvert = trired(tris(1),:);  
        centr_fir = [mean(p(indvert)) mean(q(indvert)) mean(r(indvert))];         
        distri1sq =  (centr_fir(1)-centr_i(1))^2 +(centr_fir(2)-centr_i(2))^2+...
            (centr_fir(3)-centr_i(3))^2; 
        % distri1sq = 3*(rms(centr_fir-centr_i).^2);
        
        % center of 2nd triangle
        indvert = trired(tris(2),:); 
        centr_sec = [mean(p(indvert)) mean(q(indvert)) mean(r(indvert))]; 
        distri2sq = 3*(rms(centr_sec-centr_i).^2);
        
        laplac(i,tris(1)) = 1/distri1sq ; 
        laplac(i,tris(2)) = 1/distri2sq ;
        %laplac(i,tris(3)) = 1/distri3sq ; 
        laplac(i,i) = -(1/distri1sq + 1/distri2sq);
    
    elseif length(tris) == 1 
        indvert = trired(tris(1),:);  
        centr_fir = [mean(p(indvert)) mean(q(indvert)) mean(r(indvert))];         
        distri1sq =  (centr_fir(1)-centr_i(1))^2 +(centr_fir(2)-centr_i(2))^2+...
            (centr_fir(3)-centr_i(3))^2; 
        % distri1sq = 3*(rms(centr_fir-centr_i).^2);
        
        laplac(i,tris(1)) = 1/distri1sq ; 
        %laplac(i,tris(2)) = 1/distri2sq ;
        %laplac(i,tris(3)) = 1/distri3sq ; 
        laplac(i,i) = -(1/distri1sq );
    end
end