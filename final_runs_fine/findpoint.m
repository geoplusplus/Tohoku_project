function [p,scale,angrot] = findpoint(p1,p2,somep) ; 

% [p] = findpoint(p1,p2,somep) ; 

% this finds somep in the changed refernce CS where p1 p2 on 0 ,2 on xaxis

x1 = p1(1); y1 = p1(2); 
x2 = p2(1); y2 = p2(2);

scale = 2/sqrt((x2-x1)^2+(y2-y1)^2); 
angrot = atand((y2-y1)/(x2-x1)); 

if x1>x2
    scale = -scale;
end


rotA = [cosd(angrot) sind(angrot); -sind(angrot) cosd(angrot)];
trans = scale*rotA*reshape(p1,2,1); 

p = scale*rotA*reshape(somep,2,1) - trans ; 




