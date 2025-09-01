function f=calcfs(xp,s,xB)
a = [4 2 1];
    yp = a.*xp./(sum(a.*xp));
    fs = xp(1:2)-s/(s+1)*yp(1:2)-1/(s+1)*xB(1:2);
end