function [Y_predictor, or , dmodel  , srgtSTT ]=kriging(trainX,trainY,testX,type,limit,regreeModel,relaModel,penaltyBeta,penaltyTheta)

%% 设置theta 的初始值theta0，上限ub，下限lb
lb = 10^((-1)*limit)*ones(1, size(trainX, 2)); % The lower bound of the correlation parameters
ub = 10^(limit)*ones(1, size(trainX, 2)); % The upper bound of the correlation parameters
theta0 = 1*ones(1, size(trainX, 2)); % The initial values of the correlation parameters

%% 模型
if  nargin < 7
    error('Too few input parameters')
elseif nargin ==7
    % OK UK
    [Y_predictor , or , dmodel  , srgtSTT]...
        = UK(trainX,trainY,testX,type,regreeModel,relaModel,theta0,lb,ub);    %UK
elseif  nargin > 7
    error('Too many input parameters')
end
end



