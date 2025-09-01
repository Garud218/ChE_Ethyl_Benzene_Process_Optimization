function xp =csf(xB,s,xp0)
a = [4 2 1];
flag =1; count=0;
alpha=0.2;
xp=xp0;
while flag
    fs=calcf(xp,s,xB);
    if sum(abs(fs))<1e-10
        flag=0;
    else
        Js=calcJs(xp,s,xB);
        dx=-alpha*inv(Js)*fs;
        xp(1:2)=xp(1:2)+dx;
        xp(3)=1-xp(1)-xp(2);
        count=count+1;
    end
    if count>200
        flag=0;
        disp('didnt converge in 200 iterations');
    end
end
