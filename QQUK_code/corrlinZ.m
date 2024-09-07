function  [r, dr] = corrlinZ(type,theta, d)
%CORRLIN  Linear correlation function,
%
%           n
%   r_i = prod max(0, 1 - theta_j * d_ij) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all theta_j = theta .
%
% Call:    r = corrlin(theta, d)
%          [r, dr] = corrlin(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update April 12, 2002

[n m] = size(d);  % number of differences and dimension of data
if  length(theta) == 1
  theta = repmat(theta,1,m);
% elseif  length(theta) ~= m
%   error(sprintf('Length of theta must be 1 or %d',m))
end

td=zeros(n,m);
for i=1:n
    row_z=find(type==0);
    row_x=find(type~=0);
    for j=1:m
        if ismember(j,row_x) %判断是定量因子
            td(i,j) = max(1 - abs(d(i,j)) .*(theta(j)), 0); %定量因子的高斯结构
        elseif ismember(j,row_z)  %定性因子
            if d(i,j)~=0  %定性因子的I结构
                I=1;
            else
                I=0;
            end
            td(i,j) = max(1-abs(I)*theta(j), 0);
        end
    end
end


r = prod(td, 2);

if  nargout > 1
  dr = zeros(n,m);
  for  j = 1 : m
    dr(:,j) = prod(td(:,[1:j-1 j+1:m]),2) .* (-theta(j) * sign(d(:,j)));
  end
end