function inputArgs = initializeSetup(initialCalibration)

sigmoidal = @(A,x)(A(1)./(1+exp(-(A(2).*x - A(3)))));
gaussian  = @(B,x)(B(1)*exp(-(B(2).*x - B(3)).^2));

initialCalibration.NewPolypRate = sigmoidal(initialCalibration.NewPolypRateParams,1:length(initialCalibration.NewPolypRate));
initialCalibration.EarlyProgressionRate = gaussian(initialCalibration.EarlyProgressionRateParams,1:length(initialCalibration.EarlyProgressionRate));
initialCalibration.AdvancedProgressionRate = gaussian(initialCalibration.AdvancedProgressionRate,1:length(initialCalibration.AdvancedProgressionRate));

%setting individual risks
Values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 97, 100]*5;

initialCalibration.IndividualRiskMesh = initialCalibration.IndividualRisk(Values);
% individual polyp risk
tmp(1:Values(1)) = initialCalibration.IndividualRiskMesh(1);
for x1=1:length(Values)-1
    Start = Values(x1)+1;
    Ende  = Values(x1+1);
    for x2=Start:Ende
        tmp(x2) = (initialCalibration.IndividualRiskMesh(x1) * (Ende-x2) + ...
           initialCalibration.IndividualRiskMesh(x1+1) * (x2-Start))/(Ende-Start);
    end
end
initialCalibration.IndividualRisk = tmp; clear tmp

%%

inputArgs.stats = readBenchmarkStatistics();
inputArgs.stats.mortalityYears = 5; %how many years to take into account

inputArgs.p  = 10;   % types of polyps
n  = initialCalibration.Number_patients; %%number patients

% Direct Cancer - we interpolate values (TO BE DONE BETTER!)
counter = 1;
for x1=1:19
    for x2=1:5
        inputArgs.DirectCancerRate(1, counter) = (initialCalibration.DirectCancerRate(1, x1) * (5-x2) + ...
            initialCalibration.DirectCancerRate(1, x1+1) * (x2-1))/4;
        inputArgs.DirectCancerRate(2, counter) = (initialCalibration.DirectCancerRate(2, x1) * (5-x2) + ...
            initialCalibration.DirectCancerRate(2, x1+1) * (x2-1))/4;
        counter = counter + 1;
    end
end
inputArgs.DirectCancerRate(1, counter : 150) = initialCalibration.DirectCancerRate(1, end);
inputArgs.DirectCancerRate(2, counter : 150) = initialCalibration.DirectCancerRate(2, end);
inputArgs.DirectCancerSpeed = initialCalibration.DirectCancerSpeed;

inputArgs.DirectCancerSpeed = initialCalibration.DirectCancerSpeed;

% StageVariables
inputArgs.StageVariables.Progression          = initialCalibration.Progression;
inputArgs.StageVariables.FastCancer           = initialCalibration.FastCancer;
inputArgs.StageVariables.FastCancer(6:10)     = 0;
inputArgs.StageVariables.Healing              = initialCalibration.Healing;
inputArgs.StageVariables.Symptoms             = initialCalibration.Symptoms;
inputArgs.StageVariables.Colo_Detection       = initialCalibration.Colo_Detection;
inputArgs.StageVariables.RectoSigmo_Detection = initialCalibration.RectoSigmo_Detection;
inputArgs.StageVariables.Mortality            = initialCalibration.Mortality;

inputArgs.DwellSpeed        = 'Slow';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Location                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location
inputArgs.Location.NewPolyp            = initialCalibration.Location_NewPolyp;
inputArgs.Location.DirectCa            = initialCalibration.Location_DirectCa;
inputArgs.Location.EarlyProgression    = initialCalibration.Location_EarlyProgression;
inputArgs.Location.AdvancedProgression = initialCalibration.Location_AdvancedProgression;
inputArgs.Location.CancerProgression   = initialCalibration.Location_CancerProgression;
inputArgs.Location.CancerSymptoms      = initialCalibration.Location_CancerSymptoms;
inputArgs.Location.ColoDetection       = initialCalibration.Location_ColoDetection;
inputArgs.Location.RectoSigmoDetection = initialCalibration.Location_RectoSigmoDetection;
% Location.RectoSigmoDetection = initialCalibration.Location_ColoDetection(1:10); % change later
inputArgs.Location.ColoReach           = initialCalibration.Location_ColoReach;
inputArgs.Location.RectoSigmoReach     = initialCalibration.Location_RectoSigmoReach;
% Location.RectoSigmoReach     = initialCalibration.Location_ColoReach;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   male or female                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% female 1=male, 2=female.

inputArgs.female.fraction_female               = initialCalibration.fraction_female;
inputArgs.female.new_polyp_female              = initialCalibration.new_polyp_female;
inputArgs.female.early_progression_female      = initialCalibration.early_progression_female;
inputArgs.female.advanced_progression_female   = initialCalibration.advanced_progression_female;
inputArgs.female.symptoms_female               = initialCalibration.symptoms_female;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Costs                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost

inputArgs.Cost.Colonoscopy                    = initialCalibration.Cost.Colonoscopy;
inputArgs.Cost.Colonoscopy_Polyp              = initialCalibration.Cost.Colonoscopy_Polyp;
inputArgs.Cost.Colonoscopy_Cancer             = initialCalibration.Cost.Colonoscopy_Cancer;
inputArgs.Cost.Sigmoidoscopy                  = initialCalibration.Cost.Sigmoidoscopy;
inputArgs.Cost.Sigmoidoscopy_Polyp            = initialCalibration.Cost.Sigmoidoscopy_Polyp;
inputArgs.Cost.Colonoscopy_Perforation        = initialCalibration.Cost.Colonoscopy_Perforation;
inputArgs.Cost.Colonoscopy_Serosal_burn       = initialCalibration.Cost.Colonoscopy_Serosal_burn;
inputArgs.Cost.Colonoscopy_bleed              = initialCalibration.Cost.Colonoscopy_bleed;
inputArgs.Cost.Colonoscopy_bleed_transfusion  = initialCalibration.Cost.Colonoscopy_bleed_transfusion;
inputArgs.Cost.FOBT                           = initialCalibration.Cost.FOBT;
inputArgs.Cost.I_FOBT                         = initialCalibration.Cost.I_FOBT;
inputArgs.Cost.Sept9_HighSens                 = initialCalibration.Cost.Sept9_HighSens;
inputArgs.Cost.Sept9_HighSpec                 = initialCalibration.Cost.Sept9_HighSpec;
inputArgs.Cost.other                          = initialCalibration.Cost.other;

% current treatment costs
inputArgs.CostStage.Initial(1)  = initialCalibration.Cost.Initial_I;
inputArgs.CostStage.Initial(2)  = initialCalibration.Cost.Initial_II;
inputArgs.CostStage.Initial(3)  = initialCalibration.Cost.Initial_III;
inputArgs.CostStage.Initial(4)  = initialCalibration.Cost.Initial_IV;
inputArgs.CostStage.Cont(1)     = initialCalibration.Cost.Cont_I;
inputArgs.CostStage.Cont(2)     = initialCalibration.Cost.Cont_II;
inputArgs.CostStage.Cont(3)     = initialCalibration.Cost.Cont_III;
inputArgs.CostStage.Cont(4)     = initialCalibration.Cost.Cont_IV;
inputArgs.CostStage.Final(1)    = initialCalibration.Cost.Final_I;
inputArgs.CostStage.Final(2)    = initialCalibration.Cost.Final_II;
inputArgs.CostStage.Final(3)    = initialCalibration.Cost.Final_III;
inputArgs.CostStage.Final(4)    = initialCalibration.Cost.Final_IV;
inputArgs.CostStage.Final_oc(1) = initialCalibration.Cost.Final_oc_I;
inputArgs.CostStage.Final_oc(2) = initialCalibration.Cost.Final_oc_II;
inputArgs.CostStage.Final_oc(3) = initialCalibration.Cost.Final_oc_III;
inputArgs.CostStage.Final_oc(4) = initialCalibration.Cost.Final_oc_IV;

% treatment costs in the near future
inputArgs.CostStage.FutInitial(1)  = initialCalibration.Cost.FutInitial_I;
inputArgs.CostStage.FutInitial(2)  = initialCalibration.Cost.FutInitial_II;
inputArgs.CostStage.FutInitial(3)  = initialCalibration.Cost.FutInitial_III;
inputArgs.CostStage.FutInitial(4)  = initialCalibration.Cost.FutInitial_IV;
inputArgs.CostStage.FutCont(1)     = initialCalibration.Cost.FutCont_I;
inputArgs.CostStage.FutCont(2)     = initialCalibration.Cost.FutCont_II;
inputArgs.CostStage.FutCont(3)     = initialCalibration.Cost.FutCont_III;
inputArgs.CostStage.FutCont(4)     = initialCalibration.Cost.FutCont_IV;
inputArgs.CostStage.FutFinal(1)    = initialCalibration.Cost.FutFinal_I;
inputArgs.CostStage.FutFinal(2)    = initialCalibration.Cost.FutFinal_II;
inputArgs.CostStage.FutFinal(3)    = initialCalibration.Cost.FutFinal_III;
inputArgs.CostStage.FutFinal(4)    = initialCalibration.Cost.FutFinal_IV;
inputArgs.CostStage.FutFinal_oc(1) = initialCalibration.Cost.FutFinal_oc_I;
inputArgs.CostStage.FutFinal_oc(2) = initialCalibration.Cost.FutFinal_oc_II;
inputArgs.CostStage.FutFinal_oc(3) = initialCalibration.Cost.FutFinal_oc_III;
inputArgs.CostStage.FutFinal_oc(4) = initialCalibration.Cost.FutFinal_oc_IV;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Complications                  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% risc

inputArgs.risc.Colonoscopy_RiscPerforation         = initialCalibration.Colonoscopy_RiscPerforation;
inputArgs.risc.Rectosigmo_Perforation              = initialCalibration.Rectosigmo_Perforation;
inputArgs.risc.Colonoscopy_RiscSerosaBurn          = initialCalibration.Colonoscopy_RiscSerosaBurn;
inputArgs.risc.Colonoscopy_RiscBleedingTransfusion = initialCalibration.Colonoscopy_RiscBleedingTransfusion;
inputArgs.risc.Colonoscopy_RiscBleeding            = initialCalibration.Colonoscopy_RiscBleeding;

inputArgs.risc.DeathPerforation                    = initialCalibration.DeathPerforation;
inputArgs.risc.DeathBleedingTransfusion            = initialCalibration.DeathBleedingTransfusion;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Special scenarios              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% flag, SpecialText

% we pad the strings to ensure length is always equal for the mex function
SpecialText = initialCalibration.SpecialText;
tmp = '                         ';
if length(SpecialText) >= 25
    SpecialText = SpecialText (1:25);
else tmp(1:length(SpecialText)) = SpecialText;
    SpecialText = tmp;
end

inputArgs.flag.Polyp_Surveillance  = isequal(initialCalibration.Polyp_Surveillance, 'on');
inputArgs.flag.Cancer_Surveillance = isequal(initialCalibration.Cancer_Surveillance, 'on');
inputArgs.flag.SpecialFlag         = isequal(initialCalibration.SpecialFlag, 'on');
inputArgs.flag.Screening           = isequal(initialCalibration.Screening.Mode, 'on');
inputArgs.flag.Correlation         = isequal(initialCalibration.RiskCorrelation, 'on');

% SpecialFlags
inputArgs.flag.Schoen     = false;
inputArgs.flag.Holme      = false;
inputArgs.flag.Segnan     = false;
inputArgs.flag.Atkin      = false;
inputArgs.flag.perfect    = false;
inputArgs.flag.Mock       = false;
inputArgs.flag.Kolo1      = false;
inputArgs.flag.Kolo2      = false;
inputArgs.flag.Kolo3      = false;
inputArgs.flag.Po55       = false;
inputArgs.flag.treated    = false;
inputArgs.flag.AllPolypFollowUp = false;

% poor bowel preparation
if ~contains(SpecialText, 'PBP')
    inputArgs.flag.PBP = false;
    inputArgs.PBP.Year = -1;
    inputArgs.PBP.RepeatYear = -1;
    inputArgs.PBP.Mock = -1;
else
    inputArgs.flag.PBP = true;
    tmp  = strfind(SpecialText, 'PBP');
    [tmp2, status] = str2num(SpecialText((tmp+4):(tmp+5)));
    if status
        inputArgs.PBP.Year = tmp2;
    else
        error('PBP_year could not be deciphered')
    end
    
    tmp  = strfind(SpecialText, 'WW');
    [tmp2, status] = str2num(SpecialText((tmp+2):(tmp+3))); %#ok<*ST2NM>
    if status
        inputArgs.PBP.RepeatYear = tmp2;
    else
        if isequal(SpecialText((tmp+2):(tmp+3)), 'xx')
            inputArgs.PBP.RepeatYear = -1;
        else
            error('PBP repeat year could not be deciphered')
        end
    end
    
    if  contains(SpecialText, 'Mock')
        inputArgs.PBP.Mock = 1;
    else
        inputArgs.PBP.Mock = -1;
    end
    
end


if isequal(SpecialText(1:9), 'RS-Schoen')
    inputArgs.flag.Schoen = true;
elseif isequal(SpecialText(1:8), 'RS-Holme')
    inputArgs.flag.Holme = true;
elseif isequal(SpecialText(1:9), 'RS-Segnan')
    inputArgs.flag.Segnan = true;
elseif isequal(SpecialText(1:8), 'RS-Atkin')
    inputArgs.flag.Atkin = true;
elseif isequal(SpecialText(1:7), 'perfect')
    inputArgs.flag.perfect = true;
elseif isequal(SpecialText(1:16), 'AllPolypFollowUp')
    inputArgs.flag.AllPolypFollowUp = true;
elseif isequal(SpecialText(1:5), 'Kolo1')
    inputArgs.flag.Kolo1 = true;
elseif isequal(SpecialText(1:5), 'Kolo2')
    inputArgs.flag.Kolo2 = true;
elseif isequal(SpecialText(1:5), 'Kolo3')
    inputArgs.flag.Kolo3 = true;
elseif isequal(SpecialText(1:6), 'Po+-55')
    inputArgs.flag.Po55   = true;
elseif ~isempty(regexp(SpecialText, 'treated', 'once'))
    inputArgs.flag.treated = true;
end
if ~isempty(regexp(SpecialText, 'Mock', 'once'))
    inputArgs.flag.Mock = true;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Screening Variables          %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PercentPop, Adherence, FollowUp, y-start, y-end, interval, y after colo, specificity
inputArgs.ScreeningTest(1, 1:8) = [initialCalibration.Screening.Colonoscopy(1:2), 0, initialCalibration.Screening.Colonoscopy(3:7)];
inputArgs.ScreeningTest(2, 1:8) = initialCalibration.Screening.Rectosigmoidoscopy;
inputArgs.ScreeningTest(3, 1:8) = initialCalibration.Screening.FOBT;
inputArgs.ScreeningTest(4, 1:8) = initialCalibration.Screening.I_FOBT;
inputArgs.ScreeningTest(5, 1:8) = initialCalibration.Screening.Sept9_HiSens;
inputArgs.ScreeningTest(6, 1:8) = initialCalibration.Screening.Sept9_HiSpec;
inputArgs.ScreeningTest(7, 1:8) = initialCalibration.Screening.other;

ScreeningHandles = {'Colonoscopy', 'Rectosigmoidoscopy', 'FOBT', 'I_FOBT',...
    'Sept9_HiSens', 'Sept9_HiSpec', 'other'};
% 1: colonoscopy, 2: Rectosigmoidoscopy, 3: FOBT, 4: I_FOBT
% 5: Sept9_HiSens, 6: Sept9_HiSpec, 7: other

ScreeningMatrix = zeros(1, 1000);
Start = 1;
Summe = 0;
for f=1:length(ScreeningHandles)
    Summe = Summe + initialCalibration.Screening.(ScreeningHandles{f})(1);
end
for f=1:length(ScreeningHandles)
    if initialCalibration.Screening.(ScreeningHandles{f})(1) > 0
        Ende = Start + round(initialCalibration.Screening.(ScreeningHandles{f})(1) * 1000);
        ScreeningMatrix(Start:Ende) = f;
        Start = Ende + 1;
    end
end

% P1, P2, P3, P4, P5, P6, Ca1, Ca2, Ca3, Ca4
inputArgs.Sensitivity(3,:)        = initialCalibration.Screening.FOBT_Sens;
inputArgs.Sensitivity(4,:)        = initialCalibration.Screening.I_FOBT_Sens;
inputArgs.Sensitivity(5,:)        = initialCalibration.Screening.Sept9_HiSens_Sens;
inputArgs.Sensitivity(6,:)        = initialCalibration.Screening.Sept9_HiSpec_Sens;
inputArgs.Sensitivity(7,:)        = initialCalibration.Screening.other_Sens;

% we define polyp 1-4 as early, 5-6 as advanced
inputArgs.AgeProgression      = zeros(6, 150);
inputArgs.AgeProgression(1,:) = initialCalibration.EarlyProgression    * initialCalibration.Progression(1);
inputArgs.AgeProgression(2,:) = initialCalibration.EarlyProgression    * initialCalibration.Progression(2);
inputArgs.AgeProgression(3,:) = initialCalibration.EarlyProgression    * initialCalibration.Progression(3);
inputArgs.AgeProgression(4,:) = initialCalibration.EarlyProgression    * initialCalibration.Progression(4);
inputArgs.AgeProgression(5,:) = initialCalibration.AdvancedProgression * initialCalibration.Progression(5);
inputArgs.AgeProgression(6,:) = initialCalibration.AdvancedProgression * initialCalibration.Progression(6);

inputArgs.NewPolyp              = initialCalibration.NewPolyp;              % 1:150
inputArgs.ColonoscopyLikelyhood = initialCalibration.ColonoscopyLikelyhood; % 1:150

inputArgs.IndividualRisk = zeros(1, n);
inputArgs.RiskDistribution.EarlyRisk      = initialCalibration.EarlyRisk;
inputArgs.RiskDistribution.AdvancedRisk   = initialCalibration.AdvRisk;

inputArgs.Gender         = zeros(1, n);
for f=1:n
    % we calculate an individual polyp appearance risk per patient
    inputArgs.IndividualRisk(f) = initialCalibration.IndividualRisk(round(rand*499)+1);
    % we calculate the gender of the patient. 1 = male, 2 = female.
    if rand < initialCalibration.fraction_female
        inputArgs.Gender(f) = 2;
    else
        inputArgs.Gender(f) = 1;
    end
    inputArgs.ScreeningPreference(f) = ScreeningMatrix(round(rand*999)+1);
end

% Life Table
inputArgs.LifeTable        = initialCalibration.LifeTable;
%LifeTable = zeros(size(LifeTable));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%          STAGES                   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inputArgs.StageDuration = [ 1 0 0 0;...
    0.468	    0.532       0           0;...
    0.25        0.398	    0.352       0;
    0.162       0.22        0.275	    0.343];
%     0.45	    0.55        0           0;...
%     0.12        0.135	    0.745       0;
%     0.11        0.15        0.327	    0.413];


% StageDuration = [ 1 0 0 0;...
%     0.67	    0.33        0           0;...
%     0.19        0.178	    0.632       0;
%     0.2         0.15        0.27	    0.38];
% %     0.45	    0.55        0           0;...
% %     0.12        0.135	    0.745       0;
% %     0.11        0.15        0.327	    0.413];

inputArgs.tx1 = ...
    [0.442	0.490	0.010	0.003;
    0.413	0.515	0.017	0.006;
    0.385	0.533	0.028	0.010;
    0.716	1.091	0.083	0.032;
    0.662	1.101	0.118	0.050;
    0.913	1.645	0.243	0.111;
    0.833	1.616	0.321	0.158;
    1.004	2.087	0.546	0.288;
    0.899	1.992	0.675	0.380;
    0.996	2.344	1.012	0.605;
    1.223	3.049	1.654	1.047;
    1.670	4.396	2.960	1.979;
    1.571	4.352	3.598	2.532;
    1.233	3.587	3.604	2.663;
    0.668	2.036	2.464	1.907;
    0.405	1.289	1.864	1.508;
    0.274	0.910	1.560	1.317;
    0.231	0.800	1.615	1.420;
    0.146	0.527	1.243	1.137;
    0.123	0.461	1.267	1.204;
    0.069	0.270	0.856	0.843;
    0.059	0.236	0.863	0.881;
    0.025	0.104	0.434	0.458;
    0.021	0.091	0.434	0.473;
    0.018	0.080	0.434	0.488];

% we use this matrix to conveniently assign a location to each new polyp
inputArgs.LocationMatrix = zeros(2, 1000);
Counter = 1;
% location for new polyp
for f = 1 : 13
    Ende = round(sum(inputArgs.Location.NewPolyp(1:f))/sum(inputArgs.Location.NewPolyp)*1000);
    inputArgs.LocationMatrix(1, Counter:Ende) = f;
    Counter = Ende;
end
Counter = 1;
% location for direct cancer
for f = 1 : 13
    Ende = round(sum(initialCalibration.Location_DirectCa(1:f))/sum(initialCalibration.Location_DirectCa)*1000);
    inputArgs.LocationMatrix(2, Counter:Ende) = f;
    Counter = Ende;
end

if isequal (SpecialText(1:6), 'Po+-55')
    inputArgs.LifeTable = zeros(size(inputArgs.LifeTable));
end



end

