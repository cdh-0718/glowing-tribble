function  [f, df] = regpoly2Z(type,S)
%REGPOLY2  Second order polynomial regression function
% Call:    f = regpoly2(S)
%          [f, df] = regpoly2(S)
%
% S : m*n matrix with design sites
% f =  [1 S S(:,1)*S S(:,2)S(:,2:n) ... S(:,n)^2]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update September 4, 2002
row_z=type==0;
row_x=type~=0;
X=S(:,row_x);
Z=S(:,row_z);
[m,n] = size(X);
nn = (n+1)*(n+2)/2;  % Number of columns in f  
% Compute  f
fx = [ones(m,1) X zeros(m,nn-n-1)];
j = n+1;   q = n;
for  k = 1 : n
  fx(:,j+(1:q)) = repmat(X(:,k),1,q) .* X(:,k:n);
  j = j+q;   q = q-1;
end

fz = [ones(m,1)  Z];
X_columns = size(fx, 2);
Z_columns = size(fz, 2);
result = zeros(size(fx, 1), Z_columns * X_columns);

for i = 1:Z_columns
    result(:, (i-1)*X_columns + 1 : i*X_columns) = fx.* fz(:,i);
end
f=result;


if  nargout > 1
  df = [zeros(n,1)  eye(n)  zeros(n,nn-n-1)];
  j = n+1;   q = n; 
  for  k = 1 : n
    df(k,j+(1:q)) = [2*S(1,k) S(1,k+1:n)];
    for i = 1 : n-k,  df(k+i,j+1+i) = S(1,k); end
    j = j+q;   q = q-1;
  end
end 