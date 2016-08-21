function [logfinal,final,objfn] = likelinew(model,greens,datavector,sigmasq,GPScoord,GPSdata,plot1)
% [logfinal,final] = likelinew(model,greens,datavector,sigmasq)

model = reshape(model,length(model),1); 
datavector = reshape(datavector,length(datavector),1); 

preddata = greens*model; 
error = preddata - datavector; 

final1 = -1/(2*sigmasq)*(error'*error);
final = sigmasq^(-.5)*exp(final1);
logfinal = -.5*log(sigmasq) + final1; 

objfn = error'*error; 

numdata = length(datavector);
if plot1 == 1
    predx = preddata(1:numdata/3);
    predy = preddata(numdata/3+1:2*numdata/3);
    predz = preddata(2*numdata/3+1:numdata);
    
    errx = error(1:numdata/3);
    erry = error(numdata/3+1:2*numdata/3);
    errz = error(2*numdata/3+1:numdata);
    
    figure(3);
    subplot(311); quiver(GPScoord(:,1),GPScoord(:,2),GPSdata(:,1),GPSdata(:,2))
    title('data for xy')
    axis equal; axis([-1230 1400 -700 600])
    subplot(312); quiver(GPScoord(:,1),GPScoord(:,2),predx,predy)
    title('predicted data')
    axis equal; axis([-1230 1400 -700 600])
    subplot(313); quiver(GPScoord(:,1),GPScoord(:,2),errx,erry);
    title('error')
    axis equal; axis([-1230 1400 -700 600])
    
    figure(4);
    subplot(311); scatter(GPScoord(:,1),GPScoord(:,2),10,GPSdata(:,3))
    title('data for z')
    axis equal; axis([-1230 1400 -700 600])
    subplot(312); scatter(GPScoord(:,1),GPScoord(:,2),10,predz)
    title('predicted data')
    axis equal; axis([-1230 1400 -700 600])
    subplot(313); scatter(GPScoord(:,1),GPScoord(:,2),10,errz)
    title('error')
    axis equal; axis([-1230 1400 -700 600])
end
