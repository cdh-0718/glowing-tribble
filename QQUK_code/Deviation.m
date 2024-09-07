function value = Deviation(x)
value=[];
[n,p] = size(x);

%ÖÐÐÄ»¯Æ«²îCD
mi1=zeros(n,1);mi2=zeros(n,n);mik=0;
for i=1:n
    b1 = 1;
    for j=1:p
        b1 = b1*(1+abs(x(i,j)-0.5)/2-abs(x(i,j)-0.5)^2/2);
    end
    mi1(i) = b1;
end
for i=1:n
   for k=1:n
      b2 = 1;
      for j=1:p
         b2 = b2*(1+abs(x(i,j)-0.5)/2+abs(x(k,j)-0.5)/2-abs(x(i,j)-x(k,j))/2); 
      end
      mi2(i,k) = b2;
   end
end
CD = sqrt(abs((13/12)^p-2/n*sum(mi1)+sum(sum(mi2))/n^2));
value = [value,CD];
