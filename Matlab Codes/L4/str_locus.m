function F = str_locus(a,zf,xD1,y2)
k=1./a;
x1p = xD1*k(1)/(xD1*k(1)+y2*k(2)+(1-xD1-y2)*k(3));
x2p = y2*k(2)/(xD1*k(1)+y2*k(2)+(1-xD1-y2)*k(3));
F = (y2-x2p)/(xD1-x1p)-(zf(2)-x2p)/(zf(1)-x1p);
end
