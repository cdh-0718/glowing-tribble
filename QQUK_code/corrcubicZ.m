function  [r, dr] = corrcubicZ(type,theta, d)
%CORRCUBIC  Cubic correlation function,
%
%           n
%   r_i = prod max(0, 1 - 3(theta_j*d_ij)^2 + 2(theta_j*d_ij)^3) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all theta_j = theta.
%
% Call:    r = corrcubic(theta, d)
%          [r, dr] = corrcubic(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n  matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update June 25, 2002

[n m] = size(d);  % number of differences and dimension of data
if  length(theta) == 1
  theta = repmat(theta,1,m);
% elseif  length(theta) ~= m
%   error(sprintf('Length of theta must be 1 or %d',m))
else
  theta = theta(:).';
end

td=zeros(n,m);

row_z=find(type==0);
row_x=find(type~=0);
for i=1:n
    for j=1:m
        if ismember(j,row_x) %判断是定量因子
            td(i,j) = min(abs(d(i,j)) .* theta(j), 1); %定量因子的高斯结构
        elseif ismember(j,row_z)  %定性因子
            if d(i,j)~=0  %定性因子的I结构
                I=1;
            else
                I=0;
            end
            td(i,j) = min(abs(I) .* theta(j), 1);
        end
    end
end
% td = min(abs(d) .* repmat(theta,n,1), 1);
ss = 1 - td .* (3 - 2*td.^2);
r = prod(ss, 2);

if  nargout > 1
  dr = zeros(n,m);
  for  j = 1 : m
    dd = 6*theta(j) * sign(d(:,j)) .* td(:,j) .* (td(:,j) - 1);
    dr(:,j) = prod(ss(:,[1:j-1 j+1:m]),2) .* dd;
  end
end