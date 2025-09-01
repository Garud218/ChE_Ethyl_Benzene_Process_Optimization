function f=calcf(xp,r,xD)
a = [4 2 1];
    yp = a.*xp./(sum(a.*xp));
    f = yp(1:2)-r/(r+1)*xp(1:2)-1/(r+1)*xD(1:2);
end