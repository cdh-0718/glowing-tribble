function X = spring(x)


d = 0.05+(0.2-0.05)*x(:,1);
D = 0.5 + (1.30 - 0.5)*x(:,2);

N =2*x(:,3);

for i=1:size(x,1)
y(i,1)=(N(i)+2)*D(i)*d(i)*d(i);
end
X=[d,D,x(:,3),y];