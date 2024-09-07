function X = OTL(x)

Rb1 = 50 + 100*x(:,1);
Rb2 = 25 + 45*x(:,2);
Rc1 = 1.2 + 1.3*x(:,3);
Rc2 = 0.25 + 0.95*x(:,4);
Rf = x(:,5);
Rf(Rf==1) = 0.5; Rf(Rf==2) = 1.2; Rf(Rf==3) = 2.1; Rf(Rf==4) = 2.9;
B = x(:,6);
B(B==1) = 50; B(B==2) = 100;B(B==3) = 150; B(B==4) = 200; B(B==5) = 250; B(B==6) = 300;

Vb1 = 12*Rb2./(Rb1+Rb2);

y = B.*(Vb1+0.74).*(Rc2+9)./(B.*(Rc2+9)+Rf)...
    +11.35*Rf./(B.*(Rc2+9)+Rf)...
    +(0.74*B.*Rf./Rc1).*((Rc2+9)./(B.*(Rc2+9)+Rf));
X=[Rb1,Rb2,Rc1,Rc2,x(:,5:6),y];