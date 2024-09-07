clc
clear

limit=3;

Train_n(1,:)=[100,60,40,60,80,100,80,15];
Train_n(2,:)=[100,60,40,60,80,100,80,10]*2;
Train_n(3,:)=[100,60,40,60,80,100,80,10]*3;

TF_name={'TF1','TF2','TF30','Bending','OTL','Piston','Borehole','Spring'};
for linear=1:8
    
for ni=1:3
    clear RMSE RRMSE NRMSE RMSPE
    train_n=Train_n(ni,linear);
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
%% OK,UK,UK2模型

        T1=clock;
        [Y_UK1,or_UK1,dmodel_UK1,srgtSTT_UK1]=kriging(trainX,trainY,testX,type,limit,@regpoly1Z,@corrgaussZ);
        T2=clock;
        time_UK1=etime(T2,T1);fprintf('UK1 ')

        T1=clock;
        [Y_UK2,or_UK2,dmodel_UK2,srgtSTT_UK2]=kriging(trainX,trainY,testX,type,limit,@regpoly2Z,@corrgaussZ);
        T2=clock;
        time_UK2=etime(T2,T1);fprintf('UK2 ')

        %% Y
        Y0=[testY,Y_UK1,Y_UK2];
        time(gen,:)=[time_UK1,time_UK2];

        %% 计算误差

        rrmse=zeros(size(Y0,1),size(Y0,2)-1);
        for i=1:size(Y0,1)
            for j=2:size(Y0,2)
                rrmse(i,j-1) = sqrt( mean( (Y0(i,1)-Y0(i,j)).^2)  )/std(testY, 1);%%rrmse
            end
        end
        RRMSE(gen,:)=mean(rrmse,1);
        fprintf('\n RRMSE= ');disp( RRMSE(gen,:));

    end

    %%
    trendmean(1,:)=mean(RRMSE,1);
    %% 绘制箱线图
if ni==1
    figure(linear)
    hold on 
    X_label ={'QQ_UK1','QQ_UK2'};
   set(gcf,'Units','centimeter','Position',[5 5 28 10]);
end

if ni==1
    subplot('position',[0.06,0.15, 0.26, 0.7])
    plotbox_4_3((RRMSE),X_label,'RRMSE');
    title(['n=',num2str(train_n)],'Fontsize',12,'FontWeight','normal','FontName','Times New Roman','FontWeight','bold');
elseif ni==2
    subplot('position',[0.40,0.15, 0.26, 0.7])
    plotbox_4_3((RRMSE),X_label,'RRMSE');
    title(['n=',num2str(train_n)],'Fontsize',12,'FontWeight','normal','FontName','Times New Roman','FontWeight','bold');
elseif ni==3
    subplot('position',[0.72,0.15, 0.26, 0.7])
    plotbox_4_3((RRMSE),X_label,'RRMSE');
    title(['n=',num2str(train_n)],'Fontsize',12,'FontWeight','normal','FontName','Times New Roman','FontWeight','bold');
end
hold off
end


end


T2=clock;
disp(second_change(etime(T2,T0)));


function plotbox_4_3(data,X_label,y_label,Ylim)
mycolor3 = [
    0.862745098039216,0.827450980392157,0.117647058823529;...
    0.705882352941177,0.266666666666667,0.423529411764706;...
    0.949019607843137,0.650980392156863,0.121568627450980;...
    0.956862745098039,0.572549019607843,0.474509803921569;...
    0.231372549019608,0.490196078431373,0.717647058823529;...
    0.556862745098039,0.372549019607843,0.474509803921569;...
    0.731372549019608,0.890196078431373,0.453176470588529;...
    0.556862745098039,0.272549019607843,0.674509803921569;...
    0.461372549019608,0.490196078431373,0.717647058823529;...
    0.575686274509809,0.772549019607843,0.474509803921569;...
    0.705882352941177,0.266666666666667,0.423529411764706;...
    0.949019607843137,0.650980392156863,0.121568627450980;...
    0.656862745098039,0.322549019607843,0.474509803921569;...
    0.231372549019608,0.490196078431373,0.717647058823529;...
    0.556862745098039,0.372549019607843,0.474509803921569;...
    0.956862745098039,0.572549019607843,0.474509803921569;...
    0.231372549019608,0.490196078431373,0.717647058823529;...
    0.556862745098039,0.372549019607843,0.474509803921569;...
    0.956862745098039,0.622549019607843,0.874509803921569;...
    0.231372549019608,0.490196078431373,0.717647058823529;...
    0.556862745098039,0.672549019607843,0.474509803921569;...
    0.131372549019608,0.345019607843133,0.617647058823529;
];
%坐标区域每组变量之间的标签

%% 开始绘图
%参数依次为数据矩阵、颜色设置、标记符

box_figure = boxplot(data,'Symbol','o');
%设置线宽
boxobj = findobj(gca,'Tag','Box');
for i = 1:size(data,2)
    patch(get(boxobj(i),'XData'),get(boxobj(i),'YData'),mycolor3(i,:),'FaceAlpha',0.5,...
        'LineWidth',1.5);
end
hold on;
%% 设置坐标区域的参数

set(gca,'Linewidth',1.5); %设置坐标区的线宽

% 对刻度长度与刻度显示位置调整
set(gca, 'TickDir', 'in', 'TickLength', [.008 .008]);
% 对X轴刻度与显示范围调整
set(gca,'Xlim',[0.5 size(data,2)+0.5], 'Xtick', 1:1:size(data,2),'Xticklabel',X_label,'FontName','Times New Roman','FontWeight','bold');
if nargin == 4 
set(gca,'Ylim',Ylim);
end
set(gca,'Fontsize',12,'FontWeight','bold'); % 设置坐标区字体大小
hold off
if nargin ==3
ylabel(y_label)
end
xlabel('Model');
end