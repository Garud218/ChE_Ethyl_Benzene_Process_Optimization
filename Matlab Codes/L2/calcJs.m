function Js=calcJs(xp,s,xB)
a = [4 2 1];
C = length(xB);
dfdx=zeros(C-1,C-1);
Js=zeros(C-1,C-1);
f0=calcf(xp,s,xB);
dx=1e-6;
for i=1:C-1
    idx=zeros(C,1);
    idx(i)=1;
    xplus=xp+dx*idx;
    xplus(3)=1-xplus(1)-xplus(2);
    dfdx(:,i)=(calcfs(xplus,s,xB)-f0)/dx;
end
Js=dfdx;
end