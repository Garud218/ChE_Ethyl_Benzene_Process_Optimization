function xp =crf(xD,r,xp0)
global a
flag =1; count=0;
alpha=0.2;
xp=xp0;
while flag
    f=calcf(xp,r,xD);
    if sum(abs(flag))<1e-10
        flag=0;
    else
        J=calJ(xp,r,xD);
        dx=alpha*inv(J)*f;
        xp(1:2)=xp(1:2)+dx;
        xp(3)=1-xp(1)-xp(2);
        count=count+1;
    end
    if count>200
        flag=0;
        disp('didn't converge in 200 iterations);
    end
end
