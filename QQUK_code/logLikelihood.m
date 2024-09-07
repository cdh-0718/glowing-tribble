function logL = logLikelihood(params, X, y, dim_qual, d_lv, levels)
% This function computes the loglikelihood of GP models with quantitative and qualitative variables
% params - vector of parameters, containing scale parameters phi for continuous variables 
%          and latent variables z for qualitative variables
% X - a matrix of input data, each row is a sample
% y - a vector of response data
% dim_qual - an index array of qualitative variables
% d_lv - the dimension of latent space, 1 or 2
% returns the loglikelihood value


% without qualitative variables
if nargin == 3

    n = size(X,1);
    R = computeR(params, X,X); % correlation matrix
    
    L = chol(R,'lower');
    
    one_n = ones(1,n);
    
    R_inv_y = L'\(L\y); 
    R_inv_one = L'\(L\one_n');
    
    mu = (one_n * R_inv_y) / (one_n * R_inv_one); % mean
    R_inv_y_mu = L'\(L\(y-mu));
    sigma2 = 1/n * (y-mu)' * R_inv_y_mu; % variance
    
    logL = log(det(R))+ n*log(sigma2);
    
% with qualitative variables
elseif nargin == 6

    [n,d] = size(X);
    d_qual = length(dim_qual); % number of qualitative variables
    z = params(d-d_qual+1:end); % latent variables
    
    X1 = toLatent(X,dim_qual, z, d_lv, levels); % transform to latent space
    phi = params(1:d-d_qual); % lengthscale parameters
    R = computeR([phi,zeros([1,d_lv*d_qual])], X1,X1); % correlation matrix
    
    one_n = ones(1,n);
    L = chol(R,'lower');
    
    R_inv_y = L'\(L\y); 
    R_inv_one = L'\(L\one_n');

    mu = (one_n * R_inv_y) / (one_n * R_inv_one); % mean
    R_inv_y_mu = L'\(L\(y-mu));
    
    sigma2 = 1/n * (y-mu)' * R_inv_y_mu; % variance
    logL = (log(det(R)) + n*log(sigma2));
   
end

