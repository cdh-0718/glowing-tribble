function  [r dr] = corrgaussZ(type,theta, d)
%CORRGAUSS  Gaussian correlation function,
%
%           n
%   r_i = prod exp(-theta_j * d_ij^2) ,  i = 1,...,m
%          j=1
%
% If length(theta) = 1, then the model is isotropic:
% all  theta_j = theta .
%
% Call:    r = corrgauss(theta, d)
%          [r, dr] = corrgauss(theta, d)
%
% theta :  parameters in the correlation function
% d     :  m*n matrix with differences between given data points
% r     :  correlation
% dr    :  m*n matrix with the Jacobian of r at x. It is
%          assumed that x is given implicitly by d(i,:) = x - S(i,:), 
%          where S(i,:) is the i'th design site. 

% hbn@imm.dtu.dk  
% Last update June 2, 2002
% [~,m]=size(theta);
% [n_type,m_type]=size(type);
% if size(theta,1)~= size (type,2)
%     error(sprintf('row of theta and column of type must be same'))
% end
[n,m] = size(d);  % number of differences and dimension of data

if  length(theta) == 1
  theta = repmat(theta,1,m);
elseif  length(theta) ~= m
  error(sprintf('Length of theta must be 1 or %d',m))
end


%% 多个定性因子
% T1=clock;
%     row_z=find(type==0);
%     row_x=find(type~=0);
% td=zeros(n,m);
% for i=1:n
% 
%     for j=1:m
%         if ismember(j,row_x) %判断是定量因子
%             td(i,j) = d(i,j).^2 .* (-theta(j)); %定量因子的高斯结构
%         elseif ismember(j,row_z)  %定性因子
%             if d(i,j)~=0  %定性因子的I结构
%                 I=1;
%             else
%                 I=0;
%             end
%             td(i,j) = -I*theta(j);
%         end
%     end
% end
% T2=clock;
%     time=etime(T2,T1);

% T1=clock;
row_z=type==0;
row_x=type~=0;
td = zeros(n, m);
row_x_index = false(1, m);
row_x_index(row_x) = true;
row_z_index = false(1, m);
row_z_index(row_z) = true;
d_squared = d.^2;
d_is_zero = d ~= 0;
for i = 1:n
    for j = 1:m
        if row_x_index(j) % 判断是定量因子
            td(i, j) = -d_squared(i, j) .* theta(j); % 定量因子的高斯结构
        elseif row_z_index(j)  % 定性因子
            td(i, j) = -d_is_zero(i, j) .* theta(j); % 定性因子的I结构
        end
    end
end

% T2=clock;
%     time2=etime(T2,T1);

r = exp(sum(td,2));

if  nargout > 1
  dr = repmat(-2*theta(:).',n,1) .* d .* repmat(r,1,m);
end
 