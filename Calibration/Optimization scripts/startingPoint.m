function x0 = startingPoint(opt, initialCalibration)

    x0 = [];
    for i = 1:length(opt.variablesToFit)
        tmp = initialCalibration.(opt.variablesToFit{i});
        x0 = [x0 tmp(opt.whichElementToFit{i})];
    end

end

