clc
clear


limit=3;
Train_n(2,:)=[300,120,180,160,240,300,240,30];
TF_name={'TF1','TF2','TF3','Bending','OTL','Piston','Borehole','Spring'};
for linear=1:8
    %% 参数设置
    n_starts = 200; % number of initial points of hyperparameters超参数的初始点的数量
    clear RRMSE 
    train_n=Train_n(2,linear);
    MAX=20;

    for gen=1:MAX
        test_n=1000;     %测试集
        fprintf(['训练样本n=',num2str(train_n),',','函数',num2str(linear),',第',num2str(gen),'次迭代: '])
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

%% OK,UK,UK2模型
        T1=clock;
        [Y_OK,or_OK,dmodel_OK,srgtSTT_OK]=kriging(trainX,trainY,testX,type,limit,@regpoly0Z,@corrgaussZ);
        T2=clock;
        time_OK=etime(T2,T1);fprintf('OK ')

        T1=clock;
        [Y_UK2,or_UK2,dmodel_UK2,srgtSTT_UK2]=kriging(trainX,trainY,testX,type,limit,@regpoly2Z,@corrgaussZ);
        T2=clock;
        time_UK2=etime(T2,T1);fprintf('UK2 ')

        %% Y
        Y0=[testY,Y_OK,Y_UK2];
        time(gen,:)=[time_OK,time_UK2];
        %% 计算误差
        rrmse=zeros(size(Y0,1),size(Y0,2)-1);
        data=zeros(size(Y0,1),size(Y0,2)-1);
        for i=1:size(Y0,1)
            for j=2:size(Y0,2)
                rrmse(i,j-1) = sqrt( mean( (Y0(i,1)-Y0(i,j)).^2)  )/std(testY, 1);%%rrmse
            end
        end
       RRMSE(gen,:)=mean(rrmse,1);
        fprintf('\n RRMSE= ');disp( RRMSE(gen,:));
    end

    %%
    trendmean(linear,:)=mean(RRMSE,1);
    %% 绘制箱线图
    figure(linear)
    X_label ={'QQOK','QQUK'};
    set(gcf,'Units','centimeter','Positi',[5 5 6.8 6.8]);
    ax1=axes('Position', [0.08, 0.08, 0.9, 0.85]);
    plotbox_4_4((RRMSE),X_label,TF_name{linear},'RRMSE');
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
% X_label ={'y_LV','y_MC','y_UC','y_add','PBK2','PBLK1','PBLK2'};
%% 开始绘图
%参数依次为数据矩阵、颜色设置、标记符
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
set(gca, 'Xlim', [0.5 size(data,2)+0.5], 'Xtick', 1:1:size(data,2), 'Xticklabel', X_label,'Fontsize',10,'FontWeight','bold', 'FontName', 'Times New Roman', 'FontWeight', 'bold');
% set(gca, 'Xlim', [0.5 size(data,2)+0.5], 'Xtick', 1:1:size(data,2), 'Fontsize', 12, 'FontName', 'Times New Roman', 'FontWeight', 'bold');
title(name,'Fontsize',10,'FontWeight','normal','FontName','Times New Roman','FontWeight','bold');
if nargin == 5 
set(gca,'Ylim',Ylim);
end
% set(gca,'Fontsize',8,'FontWeight','bold'); % 设置坐标区字体大小
hold off
ylabel(y_label,'Fontsize',10,'FontWeight','bold')
xlabel('模型','FontName','宋体');
end