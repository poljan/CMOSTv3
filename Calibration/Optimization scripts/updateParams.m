function initialCalibration = updateParams(x, initialCalibration, opt)

    aux = 0;
    for i = 1:length(opt.variablesToFit)
        tmp = initialCalibration.(opt.variablesToFit{i});
        tmp(opt.whichElementToFit{i}) = x(aux + (1:length(opt.whichElementToFit{i})));
        initialCalibration.(opt.variablesToFit{i}) = tmp;
        aux = aux + length(opt.whichElementToFit{i});
    end
    
end

