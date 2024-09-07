function  [f, df] = regpoly0Z(type,S)
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
[m, n] = size(S);
f = ones(m,1);
if  nargout > 1
  df = zeros(n,1);
end
