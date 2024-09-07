function X = math4(xz)
%5个定量1个定性t=3
[n,m]=size(xz);%n行，m列
x=(xz(:,1:m-1));
z1=xz(:,m);
Y=zeros(n,1);
for i=1:n
    if z1(i)==1
        %Zakharov_Function
        dim=size(x,2);
        sum_1=0;
        for H=1:dim
            sum_1=sum_1+0.5*x(i,H);
        end
        sum_2=0;
        for H=1:dim-1
            sum_2=sum_2+0.5*x(i,H)*x(i,H);
        end
        Y(i,1)=sum(x(i,:).^2)+sum_1^2+sum_1^ 4;
    elseif z1(i)==2
        %Zakharov_Function
        dim=size(x,2);
        sum_1=0;
        for j=1:dim
            sum_1=sum_1+15*(-1)^(j-1)*x(i,j);
        end
        sum_2=0;
        for j=1:dim-1
            sum_2=sum_2+25/2*x(i,j)*x(i,j+1);
        end
        sum_3=0;
        for j=1:dim
            sum_3=sum_3+10*((-1)^j)*x(i,j)^2;
        end
        Y(i,1)=(9+sum_1+sum_2+sum_3);
    elseif z1(i)==3
        %HappyCat Function
        dim=size(x,2);
        sum_1=0;
        sum_2=0;
        for j=1:dim
            sum_1=sum_1+x(i,j)^2;
            sum_2= sum_2+x(i,j);
        end
        temp1 = (abs(sum_1-dim)).^(1/4);
        temp2 = (0.5* sum_1+ sum_2)/dim;
        Y(i,1)= temp1+ temp2+0.5;
    end
end
X=[xz,Y];