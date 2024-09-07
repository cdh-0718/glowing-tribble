clc
clear

limit=3;

Train_n=[100,40,60,60,80,100,80,15];
TF_name={'TF1','TF2','TF3','Bending','OTL','Piston','Borehle','Spring'};
for linear=1:1:8
%% 参数设置
clear RRMSE 
train_n=Train_n(linear);
MAX=20;
for gen=1:MAX
    test_n=1000;     %测试集

%% 拉丁超立方训练样本
    [train ,tr_levels,dimension_dl,dimension_dx]= trainfunction(train_n,linear);% 应用示例
    X=train(:,1:dimension_dl);
    Z=train(:,dimension_dl+1:end-1);
    trainY=train(:,end);
    trainX=[X,Z];

%% 拉丁超立方测试样本

[test ,te_levels,dimension_dl,dimension_dx]= testfunction(test_n,linear);% 应用示例
if te_levels~=tr_levels
    error('The levels of dim_qual in train and test is different!')
else 
    levels=te_levels;
end

X=test(:,1:dimension_dl);
Z=test(:,dimension_dl+1:end-1);
testY=test(:,end);
testX=[X,Z];
type= [ones(1,dimension_dl),zeros(1,dimension_dx)];

dim_qual = find(type==0); % 定性因子在第几维
train=[trainX,trainY];

%% {'EXP','EXPG','GAUSS','LIN','SPHERICAL','CUBIC','SPLINE'};
T1=clock;
[Y_OK_expZ,or_OK_expZ,dmodel_OK_expZ,srgtSTT_OK_expZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@correxpZ);
T2=clock;
time_OK_expZ=etime(T2,T1);fprintf('OK_expZ ')

T1=clock;
[Y_OK_expgZ,or_OK_expgZ,dmodel_OK_expgZ,srgtSTT_OK_expgZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@correxpgZ);
T2=clock;
time_OK_expgZ=etime(T2,T1);fprintf('OK_expg ')

T1=clock;
[Y_OK_gaussZ,or_OK_gaussZ,dmodel_OK_gaussZ,srgtSTT_OK_gaussZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@corrgaussZ);
T2=clock;
time_OK_gaussZ=etime(T2,T1);fprintf('OK_gaussZ ')

T1=clock;
[Y_OK_linZ,or_OK_linZ,dmodel_OK_linZ,srgtSTT_OK_linZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@corrlinZ);
T2=clock;
time_OK_linZ=etime(T2,T1);fprintf('OK_linZ ')

T1=clock;
[Y_OK_sphericalZ,or_OK_sphericalZ,dmodel_OK_sphericalZ,srgtSTT_OK_sphericalZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@corrsphericalZ);
T2=clock;
time_OK_sphericalZ=etime(T2,T1);fprintf('OK_sphericalZ ')

T1=clock;
[Y_OK_cubicZ,or_OK_cubicZ,dmodel_OK_cubicZ,srgtSTT_OK_cubicZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@corrcubicZ);
T2=clock;
time_OK_cubicZ=etime(T2,T1);fprintf('OK_cubicZ ')

T1=clock;
[Y_OK_splineZ,or_OK_splineZ,dmodel_OK_splineZ,srgtSTT_OK_splineZ]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@corrsplineZ);
T2=clock;
time_OK_splineZ=etime(T2,T1);fprintf('OK_splineZ \n')
%% Y
Y0=[testY,Y_OK_expZ,Y_OK_expgZ,Y_OK_gaussZ,Y_OK_linZ,Y_OK_sphericalZ,Y_OK_cubicZ,Y_OK_splineZ];
time(gen,:)=[time_OK_expZ,time_OK_expgZ,time_OK_gaussZ,time_OK_linZ,time_OK_sphericalZ,time_OK_cubicZ,time_OK_splineZ];

%% 计算误差
% rmspe=zeros(size(Y0,1),size(Y0,2)-1);
rrmse=zeros(size(Y0,1),size(Y0,2)-1);
% nrmse=zeros(size(Y0,1),size(Y0,2)-1);
% rmse =zeros(size(Y0,1),size(Y0,2)-1);
data=zeros(size(Y0,1),size(Y0,2)-1);
for i=1:size(Y0,1)
    for j=2:size(Y0,2)
%         rmse (i,j-1) =sqrt( mean( (Y0(i,1)-Y0(i,j)).^2) );
        rrmse(i,j-1) = sqrt( mean( (Y0(i,1)-Y0(i,j)).^2)  )/std(testY, 1);%%rrmse
%         nrmse(i,j-1)=sqrt( mean( (Y0(i,1)-Y0(i,j)).^2)  )/(max(testY)-min(testY));
%         rmspe(i,j-1)=sqrt(sum(((Y0(i,1)-Y0(i,j))./Y0(i,1)).^2))./size(testY, 1);
    end
end


% RMSE(gen,:)=mean(rmse,1);
RRMSE(gen,:)=mean(rrmse,1);
% NRMSE(gen,:)=mean(nrmse,1);
% RMSPE(gen,:)=mean(rmspe,1);

end

%% 绘制箱线图
figure(linear)

X_label ={'QQ_EXP','QQ_EXPG','QQ_GAU','QQ_LIN','QQ_SPH','QQ_CUB','QQ_SPL'};

% RRMSE=readmatrix([num2str(linear),'result.xlsx'],'Sheet',['RRMSE',num2str(size(trainX,1))]);
set(gcf,'Units','centimeter','Position',[5 5 17 8]); 
ax1=axes('Position', [0.08, 0.08, 0.9, 0.85]);
plotbox_4_4((RRMSE),X_label,TF_name{linear},'RRMSE');
% end
end


function plotbox_4_4(data,X_label,name,y_label,Ylim)
mycolor3 = [
    0.862745098039216,0.827450980392157,0.117647058823529;...
    0.705882352941177,0.266666666666667,0.423529411764706;...
    0.949019607843137,0.650980392156863,0.121568627450980;...
    0.956862745098039,0.572549019607843,0.474509803921569;...
    0.231372549019608,0.490196078431373,0.717647058823529;...
    0.556862745098039,0.372549019607843,0.474509803921569;...
    0.731372549019608,0.890196078431373,0.453176470588529;...
];
%坐标区域每组变量之间的标签

%% 开始绘图
%参数依次为数据矩阵、颜色设置、标记符
% set(gcf,'Units','centimeter','Position',[0 0 25 12]);
box_figure = boxplot(data,'Symbol','o');
%设置线宽
boxobj = findobj(gca,'Tag','Box');
for i = 1:size(data,2)
    patch(get(boxobj(i),'XData'),get(boxobj(i),'YData'),mycolor3(i,:),'FaceAlpha',0.5,...
        'LineWidth',1);
end
hold on;
%% 设置坐标区域的参数
% xlabel(['train=',num2str(train_n*3),'dimension=',num2str(dimension)],'Fontsize',20,'FontWeight','normal');

set(gca,'Linewidth',1); %设置坐标区的线宽
% axes('Position', [0.08, 0.05, 0.91, 0.9]);
% 对刻度长度与刻度显示位置调整
set(gca, 'TickDir', 'in', 'TickLength', [.008 .008]);
% 对X轴刻度与显示范围调整
% 设置坐标轴属性
set(gca, 'Xlim', [0.5 size(data,2)+0.5], 'Xtick', 1:1:size(data,2), 'Xticklabel', X_label,'Fontsize',10.5,'FontWeight','bold', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
% set(gca, 'Xlim', [0.5 size(data,2)+0.5], 'Xtick', 1:1:size(data,2), 'Fontsize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
title(name,'Fontsize',10.5,'FontWeight','normal','FontName','Times New Roman','FontWeight','bold');
if nargin == 5 
set(gca,'Ylim',Ylim);
end
% set(gca,'Fontsize',8,'FontWeight','bold'); % 设置坐标区字体大小
hold off
ylabel(y_label,'Fontsize',10.5,'FontWeight','bold')
xlabel('Model');
end