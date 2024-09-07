%% OK
function [y1,beta,theta,s2,MSE1,RMSE1,STD1,RRMSE1,RMSPE1,MAPE1,RRMSE2,RMSPE2,MAPE2]=OK(dataset,test,flag)
[n,m]=size(dataset);
Y=(dataset(:,m));
X=(dataset(:,1:m-1));
% X= normalize(X); 
% Y = center(Y);
% [X,Y]=dsmerge(X,Y);
% dataset=[X,Y];
%模型参数设置，无特殊情况不需修改，见说明书
theta=1*ones(1,m-1); lob =1e-2*ones(1,m-1); upb = 1e2*ones(1,m-1);
%变异函数模型为高斯模型
[dmodel, perf] = dacefit(X, Y, @regpoly, @corrgauss, theta, lob, upb);
beta=dmodel.beta;
theta=dmodel.theta;
sigma2=dmodel.sigma2;
% u=F.'/R*r-fx;
% s21=sigma2*(1+u.'/(F.'/R*F)*u-r.'/R*r);
%% 计算误差


y=test(:,m);
[n,~]=size(y);
[y1,or]= predictor(test(:,1:m-1), dmodel);
[y2,s2,ME,RMSE,STD2,RRMSE]=prediction(test,dataset,theta,beta,sigma2);
Y0=[y,y1,y2];
s2=mean(or);
SST=sum((y-mean(y)).^2);
SSR=sum((y1-mean(y)).^2);
SSE=sum((y1-mean(y1)).^2);
sst=SSR+SSE;
R2=SSR/sst;
if flag==1
y1=y2;
end
MSE1=sum((y-y1).^2)/n;
RMSE1=sqrt(sum((y-y1).^2)/n);
STD1=sqrt(sum((y-mean(y)).^2)/(n-1));
RRMSE1=RMSE1/STD1;
RMSPE1=sqrt(sum(((y-y1)./y).^2))/n;
MAPE1=sum(abs(y-y1)./y)/n;

MSE2=sum((y-y2).^2)/n;
RMSE2=sqrt(sum((y-y2).^2)/n);
STD2=sqrt(sum((y-mean(y)).^2)/(n-1));
RRMSE2=RMSE2/STD2;
RMSPE2=sqrt(sum(((y-y2)./y).^2))/n;
MAPE2=sum(abs(y-y2)./(y))/n;
end