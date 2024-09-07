function X = borehole(x)

Tu = 63070 + (115600-63070)*x(:,1);
r = 100 + (50000 - 100)*x(:,2);
Hu = 990 + (1110-990)*x(:,3);
Tl = 63.1 + (116-63.1)*x(:,4);
L = 1120 + (1680-1120)*x(:,5);
Kw = 9855 + (12045 - 9855)*x(:,6);

rw = 0.05*x(:,7);
Hl = 660 + 40*x(:,8);


y = 2*pi*Tu.*(Hu-Hl)./(log(r./rw).*(1+2*L.*Tu./log(r./rw)./rw.^2./Kw+Tu./Tl));
X=[Tu,r,Hu,Tl,L,Kw,x(:,7:8),y];