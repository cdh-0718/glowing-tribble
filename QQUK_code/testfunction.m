function [test ,t,dimension_dl,dimension_dx]= testfunction(test_n,linear)
if linear==1
    dimension_dl=4;
    dimension_dx=2;
    t=[2,2];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=math1(testx);

elseif linear==2
    dimension_dl=2;
    dimension_dx=1;
    t=5;
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=math2(testx);
elseif linear==3
    dimension_dl=4;
    dimension_dx=1;
    t=3;
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=math3(testx);
elseif linear==4
    dimension_dl=2;
    dimension_dx=1;
    t=6;
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=bending(testx);

elseif linear==5
    dimension_dl=4;
    dimension_dx=2;
    t=[4,6];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=OTL(testx);
elseif linear==6
    dimension_dl=5;
    dimension_dx=2;
    t=[3,5];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=piston(testx);
elseif linear==7
    dimension_dl=6;
    dimension_dx=2;
    t=[3,4];
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=borehole(testx);
elseif linear==8
    dimension_dl=2;
    dimension_dx=1;
    t=4;
    testx=rand(test_n,dimension_dl);
    for i=1:test_n
        for k = 1:dimension_dx
            testx(i,dimension_dl+k)=ceil(rand()*t(k));
        end
    end
    test=spring(testx);

end

end
