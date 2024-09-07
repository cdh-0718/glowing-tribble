function [train ,t,dimension_dl,dimension_dx]= trainfunction(train_n,linear)
 if linear==1
    dimension_dl=4;
    dimension_dx=2;
    t=[2,2];
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=math1(x);
elseif linear==2
    dimension_dl=2;
    dimension_dx=1;
    t=5;
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=math2(x);
elseif linear==3
    dimension_dl=4;
    dimension_dx=1;
    t=3;
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=math3(x);
elseif linear==4
    dimension_dl=2;
    dimension_dx=1;
    t=6;
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=bending(x);
elseif linear==5
    dimension_dl=4;
    dimension_dx=2;
    t=[4,6];
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=OTL(x);
elseif linear==6
    dimension_dl=5;
    dimension_dx=2;
    t=[3,5];
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=piston(x);
elseif linear==7
    dimension_dl=6;
    dimension_dx=2;
    t=[3,4];
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=borehole(x); 
elseif linear==8
    dimension_dl=2;
    dimension_dx=1;
    t=4;
    train_n=ceil(train_n/prod(t));
    x=OSLHD(train_n,t,dimension_dl,dimension_dx);
    train=spring(x);


end

