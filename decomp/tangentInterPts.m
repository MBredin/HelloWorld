function [u_pAndt_pNEW,l_qAndt_qNEW] = tangentInterPts(x,fs,u_pAndt_p,l_qAndt_q,derivMethod)

u = spline(u_pAndt_p(:,1), u_pAndt_p(:,2), 1:length(x));    %interpolate upper envelope
l = spline(l_qAndt_q(:,1), l_qAndt_q(:,2), 1:length(x));    %interpolate lower envelope

xPrime = derivApprox(x(:),fs,derivMethod)';
uPrime = derivApprox(u(:),fs,derivMethod)';
lPrime = derivApprox(l(:),fs,derivMethod)';


%REESTIMATE UPPER
u_pAndt_pNEW = zeros(size(u_pAndt_p));

%HANDLE ENDS
numL = cumsum(u_pAndt_p(:,1)<0);
numL = numL(end)+1;
u_pAndt_pNEW(1:numL,:) = u_pAndt_p(1:numL,:);
numR = cumsum(u_pAndt_p(:,1)>length(x));
numR = numR(end)+1;
u_pAndt_pNEW(end-numR+1:end,:) = u_pAndt_p(end-numR+1:end,:);

%HANDEL REST INTERPLOLATION POINTS
tm = round(u_pAndt_p(:,1));
tm = tm(tm>0);
tm = tm(tm<length(x));
for k = 2:length(tm)-1

    temp = xPrime(tm(k-1):tm(k+1));  
    [~,nn] = crossing(temp,  tm(k-1):tm(k-1)+length(temp)-1,  uPrime(tm(k)) );
    [~,i] = min(abs(nn-tm(k)));
    
    m = (x(ceil(nn(i)))-x(floor(nn(i))))/(ceil(nn(i))-floor(nn(i)));
    b = x(ceil(nn(i))) - m*ceil(nn(i));
    u_pAndt_pNEW(numL-1+k,:) = [nn(i),m*nn(i)+b];
    
end


%REESTIMATE LOWER

l_qAndt_qNEW = zeros(size(l_qAndt_q));

%HANDLE ENDS
numL = cumsum(l_qAndt_q(:,1)<0);
numL = numL(end)+1;
l_qAndt_qNEW(1:numL,:) = l_qAndt_q(1:numL,:);
numR = cumsum(l_qAndt_q(:,1)>length(x));
numR = numR(end)+1;
l_qAndt_qNEW(end-numR+1:end,:) = l_qAndt_q(end-numR+1:end,:);

%HANDEL REST INTERPLOLATION POINTS
tm = round(l_qAndt_q(:,1));
tm = tm(tm>0);
tm = tm(tm<length(x));
for k = 2:length(tm)-1

    temp = xPrime(tm(k-1):tm(k+1));
    [~,nn] = crossing(temp,  tm(k-1):tm(k-1)+length(temp)-1,  lPrime(tm(k)) );
    [~,i] = min(abs(nn-tm(k)));
    
    m = (x(ceil(nn(i)))-x(floor(nn(i))))/(ceil(nn(i))-floor(nn(i)));
    b = x(ceil(nn(i))) - m*ceil(nn(i));
    l_qAndt_qNEW(numL-1+k,:) = [nn(i),m*nn(i)+b];
    
end




