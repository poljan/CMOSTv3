function optimizationResults = performOptimization(initialCalibration, benchmarks, opt)
    
    x0 = startingPoint(opt, initialCalibration)';
    
    %Fitting to data
    optsPSO = PSOSET('MAX_ITER',15,'SWARM_SIZE',35);
            
    tic
      xF = PSO('F',x0, 0*x0, 5*(x0+0.1), optsPSO, initialCalibration, benchmarks, opt);
    toc
    
    optsOptim = optimset('Display','iter','UseParallel',true,'DiffMinChange',0.01);
    optimizationResults = lsqnonlin(@F, xF,0*x0, 5*(x0+0.1),optsOptim, initialCalibration, benchmarks, opt);
    %load('res.mat')
    %tic
    %err = F(x0, initialCalibration, benchmarks, opt);
    %toc
    %sum(err.^2)
end

