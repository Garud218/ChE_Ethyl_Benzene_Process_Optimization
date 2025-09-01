function J=calcJ(xp,r,xD)
a = [4 2 1];
C = length(xD);
dfdx=zeros(C-1,C-1);
J=zeros(C-1,C-1);
f0=calcf(xp,r,xD);
dx=1e-6;
for i=1:C-1
    idx=zeros(C,1);
    idx(i)=1;
    xplus=xp+dx*idx;
    xplus(3)=1-xplus(1)-xplus(2);
    dfdx(:,i)=(calcf(xplus,r,xD)-f0)/dx;
end
J=dfdx;
end
