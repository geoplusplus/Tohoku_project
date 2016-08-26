close all;
ind = find(post==max(post)); 
ind = ind(1);

resmpl = samplestage(ind,:);

xf = linspace(1,22,22);
yf = 11:-1:1;
[xfault,yfault] = meshgrid(xf,yf);
zfault = zeros(9,22);
p = xfault(:); q = yfault(:);

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

slip1 = resmpl(11:430);
slip2 = resmpl(431:end);

figure(1); title('maxima a posteriori model')
subplot(211); patch(p(trired'),q(trired'),repmat(slip1,3,1)); axis equal
colorbar; axis off
title('maxima a posteriori model');
subplot(212); patch(p(trired'),q(trired'),repmat(slip2,3,1)); axis equal
colorbar; axis off

figure(2); title('mean model')
resmpl = mean(samplestage);
slip1 = resmpl(11:170);
slip2 = resmpl(171:end);
subplot(211); patch(p(trired'),q(trired'),repmat(slip1,3,1)); axis equal
title('mean model');colorbar;
subplot(212); patch(p(trired'),q(trired'),repmat(slip2,3,1)); axis equal
colorbar;

figure(3); title('std model')
resmpl = std(samplestage);
slip1 = resmpl(11:170);
slip2 = resmpl(171:end);
subplot(211); patch(p(trired'),q(trired'),repmat(slip1,3,1)); axis equal
title('std model');colorbar;
subplot(212); patch(p(trired'),q(trired'),repmat(slip2,3,1)); axis equal
colorbar;

% axes('Position',[.1 .7 .05 .05])
% box on
% hist(resmpl(:,171),30)