function XY = math2(X)

x = X(:,1:end-1); t = X(:,end);
y = zeros(size(X,1),1);

y(t==1) = sin(2*pi*x(t==1,2)-pi)+7*sin(2*pi*x(t==1,1)-pi).^2;
y(t==2) = sin(2*pi*x(t==2,2)-pi)+7*sin(2*pi*x(t==2,1)-pi).^2 + 12*sin(2*pi*x(t==2,2)-pi);
y(t==3) = sin(2*pi*x(t==3,2)-pi)+7*sin(2*pi*x(t==3,1)-pi).^2 + 0.5*sin(2*pi*x(t==3,2)-pi);
y(t==4) = sin(2*pi*x(t==4,2)-pi)+7*sin(2*pi*x(t==4,1)-pi).^2 + 8.0*sin(2*pi*x(t==4,2)-pi);
y(t==5) = sin(2*pi*x(t==5,2)-pi)+7*sin(2*pi*x(t==5,1)-pi).^2 + 3.5*sin(2*pi*x(t==5,2)-pi);
XY=[x,t,y];
end