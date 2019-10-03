function [err, initialCalibration] = F(x, initialCalibration,benchmarks,opt)

warning off;

initialCalibration = updateParams(x,initialCalibration, opt);

inputArgs = initializeSetup(initialCalibration);

[data.y, data.PolypsSumm, data.IncidenceCounter] = NumberCrunchingLean(inputArgs);

[~, ~, err] = Evaluation_for_Jan_Sept_2019_lean(data, initialCalibration, benchmarks, nargout > 1); % ,Step,Iter); BM

end