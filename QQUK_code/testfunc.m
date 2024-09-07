function [test ,t,dimension_dl,dimension_dx]= testfunc(test_n,linear,flag)
if linear==1
    dimension_dl=2;
    dimension_dx=1;
    t=6;
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    y=bending(testx);

elseif linear==2
    dimension_dl=4;
    dimension_dx=2;
    t=[4,6];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    y=OTL(testx);
elseif linear==3
    dimension_dl=5;
    dimension_dx=2;
    t=[3,5];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    y=piston(testx);
elseif linear==4
    dimension_dl=6;
    dimension_dx=2;
    t=[3,4];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    y=borehole(testx);
elseif linear==5
    dimension_dl=2;
    dimension_dx=1;
    t=5;
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    y=math1(testx);
end
    e=normrnd(0,0);
if flag ==1
    y=(y)+e;
else
    y=(y);
end
test=[testx,y];