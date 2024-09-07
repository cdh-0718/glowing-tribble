function [PRESSRMS,eXV,yhatXV,predvarXV,srgtSRGTXV,srgtOPTXV]=randCV(srgtOPT,kfolds,idxFolds)


    NbPoints=length(srgtOPT.T);
    if nargin==1
        kfolds=NbPoints;
        idxFolds=[1:kfolds].';
    end

    Pbkp=srgtOPT.P;
    Tbkp=srgtOPT.T;
    folds1=mod(NbPoints,kfolds);
    NbPointsPerFold1=floor(NbPoints/kfolds)+1;
    NbPointsPerFold2=floor(NbPoints/kfolds);


    srgtSRGTXV=cell(kfolds,1);
    srgtOPTXV=cell(kfolds,1);
    predvarXV=NaN(NbPoints,1);
    yhatXV=zeros(NbPoints,1);

    srgtOPT.normPar{1}=mean(Pbkp);
    srgtOPT.normPar{2}=std(Pbkp);
    srgtOPT.normPar{3}=mean(Tbkp);
    srgtOPT.normPar{4}=std(Tbkp);

    for c1=1:kfolds
        if c1<=folds1
            idx=idxFolds([((c1-1)*NbPointsPerFold1+1):c1*NbPointsPerFold1]);
        else
            idx=idxFolds([(folds1*NbPointsPerFold1+(c1-folds1-1)*NbPointsPerFold2+1)...
            :(folds1*NbPointsPerFold1+(c1-folds1-1)*NbPointsPerFold2)+NbPointsPerFold2]);
        end

        Ptest=Pbkp(idx,:);

        Ptraining=Pbkp;Ptraining(idx,:)=[];
        Ttraining=Tbkp;Ttraining(idx)=[];


        srgtOPT.P=Ptraining;
        srgtOPT.T=Ttraining;
        eval(sprintf('srgtSRGT = srgts%sFit(srgtOPT);',srgtOPT.SRGT));

        srgtSRGTXV{c1}=srgtSRGT;
        srgtOPTXV{c1}=srgtOPT;

        eval(sprintf('yhatAux = srgts%sEvaluate(Ptest, srgtSRGT);',srgtOPT.SRGT));
        yhatXV(idx,:)=yhatAux;

        switch srgtOPT.SRGT
        case {'KRG','GP'}
            eval(sprintf('predvarAux = srgts%sPredictionVariance(Ptest, srgtSRGT);',srgtOPT.SRGT));
            predvarXV(idx,:)=predvarAux;
        case 'PRS'
            eval(sprintf('predvarAux = srgts%sPredictionVariance(Ptest, Ptraining, srgtSRGT);',srgtOPT.SRGT));
            predvarXV(idx,:)=predvarAux;
        end
    end


    eXV=yhatXV-Tbkp;
    PRESSRMS=sqrt(mean(eXV.^2));

end

