function  [f, df] = regpoly1Z(type,S)
%REGPOLY1  First order polynomial regression function
%
% Call:    f = regpoly1(S)
%          [f, df] = regpoly1(S)
%
% S : m*n matrix with design sites
% f = [1  s]
% df : Jacobian at the first point (first row in S) 

% hbn@imm.dtu.dk  
% Last update April 12, 2002
row_z=type==0;
row_x=type~=0;
X=S(:,row_x);
Z=S(:,row_z);
[m,~] = size(X);
fx = [ones(m,1)  X];
fz = [ones(m,1)  Z];
X_columns = size(fx, 2);
Z_columns = size(fz, 2);
result = zeros(size(fx, 1), Z_columns * X_columns);

for i = 1:Z_columns
    result(:, (i-1)*X_columns + 1 : i*X_columns) = fx.* fz(:,i);
end
f=result;

if  nargout > 1
        df=repmat([zeros(n,1),eye(n)],1,1,m);
end