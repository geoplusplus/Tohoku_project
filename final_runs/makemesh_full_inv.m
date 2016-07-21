function [p,q,r,trired,ang,xfault,yfault,zfault] = makemesh_full_inv(surface,model,disct)


%function [p,q,r,trired2] = makemesh3(surface,polynom,strike,maxdep)
polynom(1) = model(1);
polynom(2) = model(2);
polynom(3) = 510; 
strike = 195; 
maxdep = model(3); 
amod1 = model(4);  % point solution
amod2= model(5);  % scale
bmod1 = model(6); % point solution
bmod2= model(7);   % scale

if amod1 < 2 && amod1>0 
    amod2 = 4.5* amod2; 
end

if bmod1 < 2 && bmod1>0 
    bmod2 = 4.5* bmod2; 
end

tic 
% polyfit surface trace using 3rd degree polynomial
% p = polyfit(surface(:,1),surface(:,2),3);
% f = polyval(p,surface(:,1));
% no1 = length(f);

no_s = size(surface,1); x_surf1 = surface(1,1);

xfault = []; yfault = []; zfault = [];


for i = 1:1
    x_surf = surface(i,1);
    y_surf = surface(i,2);
    z_surf = surface(i,3);
    
    if i ==1
        pol = polynom;                %[0.0813    9.8570  518.2750];
    else
        pd = x_surf - x_surf1;
        pol = polynom; pol(3) = pol(3)+pd;       %[0.0813    9.8570   518.2750+pd];
    end
    
    %%%
    %     syms z;
    %     x = pol(1)*z^2 + pol(2)*z + pol(3);
    %     pol_diff = diff(x,z,1);
    %     integrand = sqrt(1+ (pol_diff)^2);
    %     leng = (int(integrand,z));
    %     len = matlabFunction(leng);
    %     fulllen = len(z_surf) - len(-maxdep);
    %
    %     seglen = abs(fulllen/disct) ;
    %     OPTIONS = optimset('Algorithm','levenberg-marquardt','display','off');
    %     for kl = 1:disct
    %         if kl ==1
    %             func = len(z_surf) - seglen ;
    %         else
    %             func = len(zinc(kl-1)) - seglen ;
    %         end
    %     zinc(kl) = lsqnonlin(@(z) len(z) - func, -10,[],[],OPTIONS);
    %
    %     end
    
    %%%%
    pol_diff = polyder(pol);
    integrad = @(z)sqrt((pol_diff(1).*z+pol_diff(2)).^2+1);
    fulllen = integral(integrad,-maxdep,z_surf);
    len = @(z) integral(integrad,0,z) ;
    seglen = abs(fulllen/disct) ;
    OPTIONS = optimset('Algorithm','levenberg-marquardt','display','off');
    zinc = zeros(disct+1,1);
    for kl = 1:disct
        if kl ==1
            func = integral(integrad,0,z_surf) - seglen ;
        else
            func = len(zinc(kl)) - seglen ;
        end
        zinc(kl+1) = lsqnonlin(@(z) len(z) - func, -10,[],[],OPTIONS);
        
    end
    
    zinc(1) = z_surf; %zinc = zinc(1:disct+1);
    
end

parfor i = 1:no_s
    
     x_surf = surface(i,1);
    y_surf = surface(i,2);
    z_surf = surface(i,3);
    
    if i ==1
        pol = polynom;                %[0.0813    9.8570  518.2750];
    else
        pd = x_surf - x_surf1;
        pol = polynom; pol(3) = pol(3)+pd;       %[0.0813    9.8570   518.2750+pd];
    end
    
    xinc = polyval(pol,zinc); xinc(1) = x_surf;
    
    str = mean(strike(:));
    dirn = str + 90 ; % direction of the polynomial curve
    if dirn >0 && dirn <=90
        dirn = 90- dirn;
    elseif dirn > 90 && dirn <= 270
        dirn = -(dirn - 90);
    elseif dirn >270 && dirn <= 360
        dirn = 180- (dirn- 270);
    end
    
    yinc = zeros(length(xinc),1); yinc(1) = y_surf;
    for j = 2:length(xinc)
        yinc(j) = tand(dirn)*(xinc(j)-xinc(1)) + yinc(1) ;
    end
    
    xfault(:,i) = xinc; yfault(:,i) = yinc; zfault(:,i) = zinc;
    
end

p1 = [xfault(end,3) yfault(end,3)];
p2 = [xfault(end,12) yfault(end,12)];

p1ref = findpoint(p1,p2,p1); 
p2ref = findpoint(p1,p2,p2);


diffx=[]; diffy=[];
botfaultx= []; botfaulty= [];

parfor i=4:11
    top = [xfault(1,i) yfault(1,i)];
    bottom = [xfault(end,i) yfault(end,i)]; 
    
    topref = findpoint(p1,p2,top);
    bottomref = findpoint(p1,p2,bottom);
    
    Apol = (topref(2) - bottomref(2))/(topref(1)-bottomref(1));
    Bpol  = topref(2) - Apol*topref(1); 
    
    polsolve = [amod2 -(2+amod1)*amod2 (amod1*amod2*2- Apol) -Bpol];
    solx = roots(polsolve);
    
    
    % Bpol = Y2 - Apol*X2 ; 
    %syms x 
    %eqn = amod2*x*(x-2)*(x-amod1) - Apol*x -Bpol== 0; 
    %solx = solve(eqn,x);
    
    indreal = find(solx>0 &solx<2 & isreal(solx)==1);
    
    xfault_newref = double(solx(indreal));
    if isempty(indreal)
        indimag = find(imag(solx) ==0);
        xfault_newref = double(real(solx(indimag)));
    end
    xfault_newref = min(xfault_newref);
    
    yfault_newref = Apol*xfault_newref +Bpol ;
    
    fault_new = reversefindpoint(p1,p2,[xfault_newref yfault_newref]);
    xfault_new = real(fault_new(1)); 
    yfault_new = real(fault_new(2));
    
    botfaultx(i) = xfault_new; botfaulty(i) = yfault_new; 
    diffx(i) = bottom(1) - xfault_new; diffy(i) = bottom(2) - yfault_new; 
    
    
end

for i=4:11
    xfault(end,i) = botfaultx(i); yfault(end,i) = botfaulty(i);
    for j = 2: disct
        
        xfault(j,i) = xfault(j,i)- diffx(i)* (j-2)/(disct-1) ;
        yfault(j,i) = yfault(j,i)- diffy(i)* (j-2)/(disct-1) ;
    end
end


%%%%%%%%%%%%%%%%%%%%%%%
% for i=4:11
%     top = [xfault(1,i) yfault(1,i)];
%     bottom = [xfault(end,i) yfault(end,i)]; 
%     
%     topref = findpoint(p1,p2,top);
%     bottomref = findpoint(p1,p2,bottom);
%     
%     Apol = (topref(2) - bottomref(2))/(topref(1)-bottomref(1));
%     Bpol  = topref(2) - Apol*topref(1); 
%     
%     
%     % Bpol = Y2 - Apol*X2 ; 
%     syms x 
%     eqn = amod2*x*(x-2)*(x-amod1) - Apol*x -Bpol== 0; 
%     solx = solve(eqn,x);
%     
%     
%     xfault_newref = double(solx(1));
%     if isreal(xfault_newref) == 0
%         xfault_newref = double(min(real(solx)));
%     end
%     
%     yfault_newref = Apol*xfault_newref +Bpol ;
%     
%     fault_new = reversefindpoint(p1,p2,[xfault_newref yfault_newref]);
%     xfault_new = real(fault_new(1)); 
%     yfault_new = real(fault_new(2));
%     
%     
%     diffx = bottom(1) - xfault_new; diffy = bottom(2) - yfault_new; 
%     
%     xfault(end,i) = xfault_new; yfault(end,i) = yfault_new;
%     for j = 2: 15
%         
%         xfault(j,i) = xfault(j,i)- diffx* (j-2)^2/14^2 ;
%         yfault(j,i) = yfault(j,i)- diffy* (j-2)^2/14^2 ;
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%

p1 = [xfault(end,13) yfault(end,13)];
p2 = [xfault(end,20) yfault(end,20)];

p1ref = findpoint(p1,p2,p1); 
p2ref = findpoint(p1,p2,p2);


parfor i=14:19
    top = [xfault(1,i) yfault(1,i)];
    bottom = [xfault(end,i) yfault(end,i)]; 
    
    topref = findpoint(p1,p2,top);
    bottomref = findpoint(p1,p2,bottom);
    
    Apol = (topref(2) - bottomref(2))/(topref(1)-bottomref(1));
    Bpol  = topref(2) - Apol*topref(1); 
    
    
    % Bpol = Y2 - Apol*X2 ; 
    polsolve = [bmod2 -(2+bmod1)*bmod2 (bmod1*bmod2*2- Apol) -Bpol];
    solx = roots(polsolve);
    
%     syms x 
%     eqn = bmod2*x*(x-2)*(x-bmod1) - Apol*x -Bpol == 0; 
%     solx = solve(eqn,x);
    

    indreal = find(solx>0 &solx<2 & isreal(solx)==1);
    
    xfault_newref = double(solx(indreal));
    if isempty(indreal)
        indimag = find(imag(solx) ==0);
        xfault_newref = double(real(solx(indimag)));
    end
    xfault_newref = min(xfault_newref);
    
%     xfault_newref = double(solx(1));
%      if isreal(xfault_newref) == 0
%          indimag = find(imag(solx) ==0);
%         xfault_newref = double(real(solx(indimag)));
%         %xfault_newref = double(min(real(solx)));
%     end
    yfault_newref = Apol*xfault_newref +Bpol ;
    
    fault_new = reversefindpoint(p1,p2,[xfault_newref yfault_newref]);
    xfault_new = real(fault_new(1)); 
    yfault_new = real(fault_new(2));
    
    
    diffx(i) = bottom(1) - xfault_new; diffy(i) = bottom(2) - yfault_new; 
    botfaultx(i)= xfault_new; botfaulty(i) = yfault_new;
    
%     xfault(end,i) = xfault_new; yfault(end,i) = yfault_new;
%     for j = 2: 15
%         
%         xfault(j,i) = xfault(j,i)- diffx* (j-2)^2/14^2 ;
%         yfault(j,i) = yfault(j,i)- diffy* (j-2)^2/14^2 ;
%     end
end
for i = 14:19
    xfault(end,i) = botfaultx(i); yfault(end,i) = botfaulty(i);
    for j = 2: disct
        
        xfault(j,i) = xfault(j,i)- diffx(i)* (j-2)/(disct-1) ;
        yfault(j,i) = yfault(j,i)- diffy(i)* (j-2)/(disct-1) ;
    end
end

% for i=14:19
%     top = [xfault(1,i) yfault(1,i)];
%     bottom = [xfault(end,i) yfault(end,i)]; 
%     
%     topref = findpoint(p1,p2,top);
%     bottomref = findpoint(p1,p2,bottom);
%     
%     Apol = (topref(2) - bottomref(2))/(topref(1)-bottomref(1));
%     Bpol  = topref(2) - Apol*topref(1); 
%     
%     
%     % Bpol = Y2 - Apol*X2 ; 
%     
%     syms x 
%     eqn = bmod2*x*(x-2)*(x-bmod1) - Apol*x -Bpol == 0; 
%     solx = solve(eqn,x);
%     
%     xfault_newref = double(solx(1));
%      if isreal(xfault_newref) == 0
%         xfault_newref = double(min(real(solx)));
%     end
%     yfault_newref = Apol*xfault_newref +Bpol ;
%     
%     fault_new = reversefindpoint(p1,p2,[xfault_newref yfault_newref]);
%     xfault_new = real(fault_new(1)); 
%     yfault_new = real(fault_new(2));
%     
%     
%     diffx = bottom(1) - xfault_new; diffy = bottom(2) - yfault_new; 
%     
%     xfault(end,i) = xfault_new; yfault(end,i) = yfault_new;
%     for j = 2: 15
%         
%         xfault(j,i) = xfault(j,i)- diffx* (j-2)^2/14^2 ;
%         yfault(j,i) = yfault(j,i)- diffy* (j-2)^2/14^2 ;
%     end
% end

zdiff = abs(mean(zfault(end,:)) - mean(zfault(end-1,:)));
zdiff1 = abs(mean(zfault(end-1,:)) - mean(zfault(end-2,:)));
xdiff = abs(mean(xfault(end,:)) - mean(xfault(end-1,:)));
xdiff1 = abs(mean(xfault(end-1,:)) - mean(xfault(end-2,:)));

theta = atand(zdiff/xdiff);
theta1 = atand(zdiff1/xdiff1);

diffall = sqrt(diffx.^2+diffy.^2);

if theta> 40 || max(diffall) > 50 || theta1 > 35
    ang = 0;
else ang =1;
end

disct = disct; col = no_s;
%  figure; surf(xfault,yfault,zfault); axis equal
xfault = xfault(:,1:2:end); yfault = yfault(:,1:2:end); zfault = zfault(:,1:2:end);

p = xfault(:); q= yfault(:); r= zfault(:);

rowp = size(xfault,1); 
colp = size(xfault,2);
numpatch = 2*(rowp-1)*(colp-1);

trired = []; 
for i = 1: colp-1
    for j = 1: rowp-1
        trired = [trired; (i-1)*rowp+j (i-1)*rowp+j+1 i*rowp+j+1];
        trired = [trired; (i-1)*rowp+j i*rowp+j i*rowp+j+1];
    end
end


%trimesh(trired,p,q,r)






