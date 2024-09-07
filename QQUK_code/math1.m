function X = math1(xz)
%1个定量1个定性t=4

[n,m]=size(xz);%n行，m列
x=(xz(:,1:m-2));
z1=xz(:,m-1:m);
Y=zeros(n,1);
x(:,1)=1+x(:,1);
x(:,2)=1+x(:,2);
x(:,3)=2*pi*x(:,3);
x(:,4)=2*pi*x(:,4);
for i=1:n
    if z1(i,1)==1&&z1(i,2)==1
        Y(i,1) = x(i,1)+sin(x(i,2))+x(i,3)+exp(x(i,4));
    elseif z1(i,1)==1&&z1(i,2)==2
        Y(i,1)= x(i,1)-cos(x(i,2))+x(i,3)-exp(x(i,4));
    elseif z1(i,1)==2&&z1(i,2)==1
        Y(i,1)= x(i,1)*sin(x(i,2))+x(i,3)*exp(x(i,4));
    elseif z1(i,1)==2&&z1(i,2)==2
        Y(i,1)= cos(x(i,1))/x(i,2)+exp(x(i,3))/x(i,4);
    elseif z1(i,1)==3&&z1(i,2)==1
        Y(i,1)= x(i,1)*sin(x(i,2))+x(i,3)*exp(x(i,4));
    elseif z1(i,1)==3&&z1(i,2)==2
        Y(i,1)= x(i,1)*cos(x(i,2))+x(i,3)*exp(x(i,4));
    end
end
X=[x,z1,Y];
