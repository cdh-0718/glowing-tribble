% clear all
% m=50;
% t=3;
% q=4;

function Data=OSLHD(m,t1,dimension_dl,dimension_dx)
%m是每层的数量
%t是层数
%dimension_dl是定量因子个数，dimension_dx是定性因子个数
% t=t1^dimension_dx;
t=1;
for dim=1:length(t1)
    t=t*t1(dim);
end
if dim~=dimension_dx
    error('t and dimension_dx must have the same number of rows'), end
% m=ceil(m/t);

n=m*t;
X=[];
for j=1:t
    x{j}=[];
    for i=1:dimension_dl
        Data=randperm(m).';
        x{j}=[x{j},Data];
    end

    X=[X;x{j}];
    X1=X;
end
flag=zeros(n,dimension_dl);
for j=1:dimension_dl
    for l=1:m
        g=1;
        Nu=randperm(t)+t*(l-1);
        for i=1:n
            if X(i,j)==l&&flag(i,j)==0
                X(i,j)=Nu(g);
                g=g+1;
                flag(i,j)=l;
                F(i,j)=l;
            end
        end
    end
end
Markov_length =500;
Energy_current = inf;                        % 当前解的能量
Energy_best = inf;
beta=0.1;
res=1;
temperature=500;
while temperature > res
    Energy1=Energy_best;                             % 用于控制循环的结束条件

    for h = 1: Markov_length
        A=rand();
        % 产生新解(对当前解添加扰动)
        if A<=0.5                   % 切片中某列两个元素交换
            z=randperm(t);
            for j=1:t
                x{j}=X((j-1)*m+1:m*j,:);
            end
            randa=ceil(m*rand());
            randb=ceil(m*rand());
            b=ceil(dimension_dl*rand());
            x{z(1)}([randa,randb],b)=x{z(1)}([randb,randa],b);
            Fi=(z(1)-1)*m+randa;
            Fj=(z(1)-1)*m+randb;
            F([Fi,Fj],b)=F([Fj,Fi],b);
            X=[];
            for j=1:t
                X=[X;x{j}];
            end
            %  clear z randa randb b j l Nu flag
        else

            j=ceil(dimension_dl*rand());
            l=ceil(m*rand());
            g=1;
            Nu=randperm(t)+t*(l-1);
            flag=zeros(n,dimension_dl);
            for i=1:n
                if F(i,j)==l && flag(i,j) ==0
                    X(i,j)=Nu(g);
                    g=g+1;
                    flag(i,j)=l;
                    F(i,j)=l;
                end
            end
        end
        % 回到起始点，加上首尾两个城市的距离
        route_new=X;
        X_Energy =  Deviation(X);
        for i=1:t
            route{i}=route_new((i-1)*m+1:i*m,:);
            Energy(i)= Deviation(route{i});
        end
        x_Energy=mean(Energy,2);
        Energy_new=mean(X_Energy+x_Energy);


        % 按照Metroplis准则接收新解
        if Energy_new<Energy_current
            % 更新局部最优
            Energy_current = Energy_new;
            route_current = route_new;
            % 更新全局最优
            if Energy_new<Energy_best
                Energy_best=Energy_new;
                route_best = route_new;
            end
        else
            if rand<exp(-(Energy_new-Energy_current)/temperature)
                Energy_current = Energy_new;
                route_current = route_new;
            else
                route_new = route_current; % 否则路线不更新，保存更改之前的路线
            end
        end
        %             fprintf(' 中心化偏差CD = %5.4f , temperature = %5.4f / %d .\n',Energy_best,temperature,res);
    end

    %     if Energy1==Energy_best&&Energy_current==Energy_best
    %         break
    %     else
    %         temperature = temperature*ratio; %降温过程
    temperature = temperature/(1+beta*temperature); %降温过程


end

X(:,1:dimension_dl)=route_best(:,1:dimension_dl)./(m*t);
Data_dl=X;
clear X Data
%% 定性因子的赋值
Z= fullfact(t1);
Data_dx = repmat(Z, m, 1);

columnorder = 1:dimension_dx;  % 创建从1到dimension_dx的行向量
columnorder = flip(columnorder);  % 反转行向量得到从dimension_dx到1的排列
Data_dx = sortrows(Data_dx, columnorder);% dimension_dx从最后一列往前依次升序

Data=[Data_dl,Data_dx];
end
