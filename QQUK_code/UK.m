function[Y,Y_or,dmodel,perf]=UK(trainX,trainY,testX,type,regreeModel,relaModel,theta0,lb,ub)

[dmodel, perf]=dacefit_UK(trainX, trainY, type,regreeModel, relaModel, theta0, lb, ub);
[Y,Y_or]= predictor(testX,dmodel);
end