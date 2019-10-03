% This script can be used to calibrate CMOST for given benchmarks

clear variables; %clearing variables present in the workspace

addpath('Additional functions')
addpath('Optimization scripts')
addpath('../Additional procedures')
addpath('../Core procedures');

%first load benchmarks for the calibration
benchmarks = loadBenchmarks([]);

%load initial calibration
[initialCalibration, path, file] = readCalibration([]);

%%
%adding parametrizations of functions
initialCalibration.NewPolypRateParams = [0.0267, 0.4925, 6.5697]; %a0, a1, a2
initialCalibration.EarlyProgressionRateParams = [0.1109, 0.0668, 0.422]; %b0, b1, b2
initialCalibration.AdvancedProgressionRateParams = [0.0158, 0.0850, 2.1862]; %c0, c1, c2
initialCalibration.IndividualRiskMesh = [0.0701, 0.1402, 0.2104, 0.4207,  0.7713, 1.4725,  2.6996,  3.6711,  4.1366,  4.5044, 23.8956, 52.9823];

initialCalibration.Number_patients = 1000000; %number of patients in each run, the bigger the better

% %define which parameters are to be fiddled with
opt.variablesToFit = {'new_polyp_female','Progression','early_progression_female','advanced_progression_female',...
                    'NewPolypRateParams','IndividualRiskMesh',...
				  'EarlyProgressionRateParams',...
				  'Location_EarlyProgression','Location_AdvancedProgression',...
				  'FastCancer','AdvancedProgressionRateParams'};
              
% parameters can be vectors so I need to specify which elements of each
% parameter to fit, order as in the definition above
opt.whichElementToFit = {[1],[1:3 5],[1],[1],...
    [1:3],[1:12],...
    [1:3],...
    [13],[13],...
    [1:5],[1:3]}; %this indicate which elements of given variables to fit should be tweaked

%% 
%perform optimization
optimizationResults = performOptimization(initialCalibration, benchmarks, opt);

%% plot (with comparizon with initial calibration) nad save the results
[~, optimizedCalibration] = F(optimizationResults, initialCalibration, benchmarks, opt);

calibration.Variables = optimizedCalibration;
save([path 'optimized' datestr(datetime(),'MM_dd_yy_hh_mm_ss') '_' file],'calibration');

