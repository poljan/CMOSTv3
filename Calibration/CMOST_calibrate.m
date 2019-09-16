%% This script can be used to calibrate CMOST for given benchmarks

clear variables; %clearing variables present in the workspace

addpath('Additional functions')
addpath('../Additional procedures')
addpath('../Core procedures');

%first load benchmarks for the calibration
benchmarks = loadBenchmarks([]);

%load initial calibration
initialCalibration = readCalibration([]);

inputArgs = initializeSetup(initialCalibration);

%%
%define which parameters are to be fiddled with
%Progression variable is a vector and we can fit each element separately
opt.variablesToFit = {'new_polyp_female','Progression','early_progression_female','advanced_progression_female'
                    'NewPolypRate','IndividualRisk',... %those were for step 1
				  'EarlyProgressionRate',
                  ,... %those were for step 2
				  ,'Location_EarlyProgression','Location_AdvancedProgression',...
				  'FastCancer','AdvancedProgressionRate'};
              
              %I need to parametrize those functions
opt.whichElementToFit = {[1], ...
                        [1:3 5],[1:20],[1:500],[1:20],...
                        [1],[1]}; %this indicate which elements of given variables to fit should be tweaked
opt.lowerVariableBounds = {};
opt.upperVariableBounds = {};

%perform optimization
%optimizationResults = performOptimization(initialCalibration, benchmarks, opt);

%save and plot the results
