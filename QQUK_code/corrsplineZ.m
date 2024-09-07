function  [r, dr] = corrsplineZ(type,theta, d)
%CORRSPLINE  Cubic spline correlation function,
%
%           n
%   xi = prod S(theta_j*d_ij) ,  i = 1,...,m
%          j=1
%
% with
%           1 - 15xi^2 + 30xi^3   for   0 <= xi <= 0.2
%   S(xi) =  1.25(1 - xi)^3       for  0.2 < xi < 1
%           0                   for    xi >= 1
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta.
%
% Call:    r = corrspline(theta, d)
%          [r, dr] = corrspline(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update May 30, 2002

[n m] = size(d);  % number of differences and dimension of data
if  length(theta) == 1
  theta = repmat(theta,1,m);
% elseif  length(theta) ~= m
%   error(sprintf('Length of theta must be 1 or %d',m))
else
  theta = theta(:).';
end
mn = n*m;   ss = zeros(mn,1);

td=zeros(n,m);
row_z=find(type==0);
row_x=find(type~=1);
for i=1:n

    for j=1:m
        if ismember(j,row_x) %判断是定量因子
            td(i,j) = abs(d(i,j)) .* theta(j); %定量因子的高斯结构
        elseif ismember(j,row_z)  %定性因子
            if d(i,j)~=0  %定性因子的I结构
                I=1;
            else
                I=0;
            end
            td(i,j) = abs(I) .* theta(j);
        end
    end
end
xi = reshape(td, mn,1);
% xi = reshape(abs(d) .* repmat(theta,n,1), mn,1);
% Contributions to first and second part of spline
i1 = find(xi <= 0.2);
i2 = find(0.2 < xi & xi < 1);
if  ~isempty(i1)
  ss(i1) = 1 - xi(i1).^2 .* (15  - 30*xi(i1));
end
if  ~isempty(i2)
  ss(i2) = 1.25 * (1 - xi(i2)).^3;
end
% Values of correlation
ss = reshape(ss,n,m);
r = prod(ss, 2);

if  nargout > 1  % get Jacobian
  u = reshape(sign(d) .* repmat(theta,n,1), mn,1);
  dr = zeros(mn,1);
  if  ~isempty(i1)
    dr(i1) = u(i1) .* ( (90*xi(i1) - 30) .* xi(i1) );
  end
  if  ~isempty(i2)
    dr(i2) = -3.75 * u(i2) .* (1 - xi(i2)).^2;
  end
  ii = 1 : n;
  for  j = 1 : m
    sj = ss(:,j);  ss(:,j) = dr(ii);
    dr(ii) = prod(ss,2);
    ss(:,j) = sj;   ii = ii + n;
  end
  dr = reshape(dr,n,m);
end  % Jacobian