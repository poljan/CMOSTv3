function [currentYear, PolypsSumm, IncidenceCounter] = NumberCrunchingLean(p, StageVariables, Location, Cost, CostStage, risc,...
    flag, SpecialText, female, Sensitivity, ScreeningTest, ScreeningPreference, AgeProgression,...
    NewPolyp, ColonoscopyLikelyhood, IndividualRisk, RiskDistribution, Gender, LifeTable,...
    LocationMatrix, StageDuration, tx1, DirectCancerRate, DirectCancerSpeed,DwellSpeed, PBP, stats)

inputArgs = p;
names = fieldnames(inputArgs);
for i=1:length(names)
    eval([names{i} '= inputArgs.' names{i} ';' ]);
end

simSettings.yearsToSimulate = 100; %numbers of years to simulate
simSettings.numPeriods = 4; %number of periods in each simulation year

NAlive = uint32(length(Gender)); %number of individuals
%disp(['Simulated population size: ' int2str(NAlive)]);
SubjectIDs = uint32(1:NAlive)'; %this will hold actual indices to the output vector
NumIndividualsToSimulate = NAlive;

BlockFromDeath = false(NAlive,1);

%substitute of Included
SubjectsIDWouldBeAlive = uint32([]);
GenderWouldBeAlive     = uint16([]);

%precomputing life table by taking into account number of simulated periods
%within a year
LifeTable = 1-(1-LifeTable').^(1/double(simSettings.numPeriods));

%recasting other input variables
Gender = uint16(Gender)';
IndividualRisk = single(IndividualRisk');
ScreeningPreference = uint8(ScreeningPreference');
NewPolyp = single(NewPolyp');
RiskDistribution.EarlyRisk = single(RiskDistribution.EarlyRisk');
RiskDistribution.AdvancedRisk = single(RiskDistribution.AdvancedRisk');

Counter = 1;
for f = 1 : 13
    Ende = round(sum(Location.NewPolyp(1:f))/sum(Location.NewPolyp)*1000);
    LocationMatrix(Counter:Ende) = f;
    Counter = Ende;
end
LocationMatrix = uint8(LocationMatrix'); %again we can use uint8 for that

flds = fields(StageVariables);
for i = 1:length(flds) %transpose each field
    StageVariables.(flds{i}) = StageVariables.(flds{i})';
end
StageVariables.FastCancerRatios = StageVariables.FastCancer;
for i = 2:10
    StageVariables.FastCancerRatios(i) = StageVariables.FastCancer(i)/StageVariables.FastCancer(i-1);
end

flds = fields(CostStage);
for i = 1:length(flds) %transpose each field
    CostStage.(flds{i}) = CostStage.(flds{i})';
end

%recasting control variables
yearsToSimulate = uint16(simSettings.yearsToSimulate);
numPeriods = uint16(simSettings.numPeriods);

%preparing additional variables

%matrix for fast indexing
%here again we don't need to have double precision
GenderProgression = single(ones(10, 2));
GenderProgression(1:4, 2) = single(female.early_progression_female);
GenderProgression(5:6, 2) = single(female.advanced_progression_female);
GenderProgressionRows = uint16(size(GenderProgression,1));
GenderProgressionRatios = GenderProgression;
for i = 2:size(GenderProgression,1)
    GenderProgressionRatios(i,:) =  GenderProgression(i,:)/GenderProgression(i-1,:);
end

% matrix for fast indexing
LocationProgression = zeros(10, 13); % 10 polyps x 13 locations
LocationProgression(1:5, :) = repmat(Location.EarlyProgression,5,1);
LocationProgression(6, :) = Location.AdvancedProgression;
LocationProgression(7:10, :) = repmat(Location.CancerProgression,4,1);
LocationProgressionRatios = LocationProgression;
for i = 2:10
    LocationProgressionRatios(i,:) = LocationProgression(i,:)./LocationProgression(i-1,:);
end

LocationProgression = single(LocationProgression'); %again we don't need double precision
LocationProgressionRatios = single(LocationProgressionRatios'); %again we don't need double precision

LocationProgressionRows = uint16(size(LocationProgression,1)); %number of rows

% Cancer progression
%
StageDurationC = cumsum(StageDuration,2);
StageDurationC(StageDurationC == 1) = Inf;

%All of the matrices responsible for generating multinomail random
%variables are replaced by the Inverse CDF approach, i.e. by appropriate
%griddedInterpolants
%this saves memory especially in the case of Mortality time generation

%the beloow lines are depreciated, now in the main loop there is gender and
%age dependent stage random number generator
% Stage 1: 5%, Stage 2: 35%, Stage 3: 40%, Stage 4: 20% - THIS IS NOT TRUE
%StageProbability = [0 cumsum([150 356 279 215]/1000)];
%StageRandomGenerator = griddedInterpolant(StageProbability,7:11,'previous');

%defining sojourn time random generator
SojournCDF = [zeros(1,4); cumsum(bsxfun(@rdivide,tx1,sum(tx1)))]; %cumulative distribution function
meshSojournCDF = unique(SojournCDF);
SojournInvCDF = zeros(length(meshSojournCDF),4);
for i = 1:4
    SojournInvCDF(:,i) = interp1(SojournCDF(:,i)', 1:26, meshSojournCDF);
end
F = griddedInterpolant({7:10, meshSojournCDF},SojournInvCDF');
SojournRandomGenerator = @(Stage, R)(uint16(floor(F(Stage,R))));

% reach of colonoscopy
F = find(Location.ColoReach == 1, 1,'first');
GI = griddedInterpolant([0 Location.ColoReach(1:F)],1:(F+1),'previous');
ColoReachRandomGenerator = @(R)(uint8(GI(R)));
Location.ColoDetection = Location.ColoDetection'; %transposing

%output variables
%arrays with polyps will be dynamically expanded and contracted
Polyp.Polyps            = uint16([]);   % we initialize the polyps cells
Polyp.PolypYear         = uint16([]);   % we initialize the field for keeping track of years
Polyp.PolypLocation     = uint16([]);   %with uint8 we have up to 255 locations
Polyp.AdvProgression    = uint16([]); %up to 65535 with uint16
Polyp.EarlyProgression  = uint16([]); %up to 65535 with uint16
Polyp.SubjectID         = uint32([]);
Polyp.ProgressionProb   = single([]);
Polyp.DirectProgressionProb = single([]);
Polyp.OwnersGender      = uint16([]);

Ca.Cancer               = uint8([]);
Ca.CancerYear           = uint16([]);
Ca.CancerLocation       = uint8([]);
Ca.DwellTime            = uint16([]);
Ca.SympTime             = uint16([]);
Ca.SympStage            = uint16([]);
Ca.TimeStage_I          = uint16([]);
Ca.TimeStage_II         = uint16([]);
Ca.TimeStage_III        = uint16([]);
Ca.SubjectID            = uint32([]);
Ca.OwnersGender         = uint16([]);

Detected.Cancer         = uint8([]);
Detected.CancerYear     = uint16([]);
Detected.CancerLocation = uint8([]);
Detected.MortTime       = uint16([]);
Detected.SubjectID      = uint32([]);

%tmp = uint16(zeros(NAlive,yearsToSimulate));
%NumCancer          = tmp;
%MaxCancer          = tmp;

%DeathYear  = uint16(zeros(NAlive,1));

% TumorRecord.Stage         = uint8([]);
% TumorRecord.Location      = uint8([]);
% TumorRecord.DwellTime     = uint16([]);
% TumorRecord.Sojourn       = uint8([]);
% TumorRecord.Gender        = uint8([]);
% TumorRecord.Detection     = uint8([]);
% TumorRecord.SubjectID     = uint32([]);
% TumorRecord.Year          = uint16([]);
% TumorRecord.Time          = uint16([]);

Last.Colonoscopy  = zeros(1,NAlive);
Last.Polyp        = ones(1,NAlive) *-100;
Last.AdvPolyp     = ones(1,NAlive) *-100;
Last.Cancer       = ones(1,NAlive) *-100;
Last.ScreenTest   = zeros(1,NAlive);

PolypsSumm.NumPolyps = zeros(6,yearsToSimulate);
PolypsSumm.FracPolyps = zeros(6,yearsToSimulate);
PolypsSumm.NumPolypsMale = zeros(6,yearsToSimulate);
PolypsSumm.FracPolypsMale = zeros(6,yearsToSimulate);
PolypsSumm.NumPolypsFemale = zeros(6,yearsToSimulate);
PolypsSumm.FracPolypsFemale = zeros(6,yearsToSimulate);

PolypsSumm.ActNumPolyps = zeros(5,yearsToSimulate);

IncidenceCounter = zeros(6, yearsToSimulate); %overall, male, female (first 3 rows is num cancers, next 3 rows are the number of indiviudals)

%% MAIN SIMULATION LOOP
stepCounter = uint16(0); %
ni = uint16(stats.mortalityYears);
while and(NAlive > 0 || ~isempty(GenderWouldBeAlive), stepCounter < yearsToSimulate*numPeriods)
    currentYear = idivide(stepCounter,numPeriods)+1; %need to use idivide, because of interger division rounding
    currentYear = uint8(currentYear); %just to be on a safe side when indexing matrices
    stepCounter = stepCounter + 1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  some precalculation steps, done only once a year %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(stepCounter,numPeriods) == 1 %if we start a new year
        
        PolypRate = IndividualRisk .* NewPolyp(currentYear);
        % the gender specific risk
        PolypRate(Gender==2) = PolypRate(Gender==2) * female.new_polyp_female;
        
        %% mortality time generator preparations
        ranges = uint8([51 60 70 80]);
        whichCurve = sum(currentYear >= ranges);
        
        %I need to select survival curves by gender and for current age
        %group
        % SurvivalTmp = Survival(:,whichCurve*4 + (1:4))';
        nEl = 11; %number of data points for survival, there is data for 0,12,..,120 months
        nAgeGroups = 5; %number of available age groups
        nStages = 4;
        
        %gender -> stage -> age
        femaleShift = nEl*nAgeGroups*nStages;
        stageShift = nEl*nAgeGroups;
        indxS = ((1:nEl) + whichCurve*nEl)';
        SurvivalTmpMale   = stats.osByGenderAgeStage.surv([indxS indxS+stageShift indxS+2*stageShift indxS+3*stageShift])';
        SurvivalTmpFemale = stats.osByGenderAgeStage.surv([indxS indxS+stageShift indxS+2*stageShift indxS+3*stageShift]+femaleShift)';
        
        %%
        
        % we create a smooth curve (we take only first 6 points)
        L = double(ni*numPeriods+1);
        MortalityCDFmale = zeros(L+1,4); %prepare matrix for CDF
        MortalityCDFfemale = zeros(L+1,4); %prepare matrix for CDF
        for f = 1:4
            MortalityCDFmale(:,f) = [1-interp1(linspace(0,1,6), SurvivalTmpMale(f,1:(ni+1)),linspace(0,1,L)) 1];
            MortalityCDFfemale(:,f) = [1-interp1(linspace(0,1,6), SurvivalTmpFemale(f,1:(ni+1)),linspace(0,1,L)) 1];
        end
        
        meshMortalityCDFmale = unique(MortalityCDFmale);
        meshMortalityCDFfemale = unique(MortalityCDFfemale);
        MortalityInvCDFmale = zeros(length(meshMortalityCDFmale),4);
        MortalityInvCDFfemale = zeros(length(meshMortalityCDFfemale),4);
        for i = 1:4
            MortalityInvCDFmale(:,i) = interp1(MortalityCDFmale(:,i)', 1:(L+1), meshMortalityCDFmale);
            MortalityInvCDFfemale(:,i) = interp1(MortalityCDFfemale(:,i)', 1:(L+1), meshMortalityCDFfemale);
        end
        FmortMale = griddedInterpolant({7:10, meshMortalityCDFmale},MortalityInvCDFmale');
        FmortFemale = griddedInterpolant({7:10, meshMortalityCDFfemale},MortalityInvCDFfemale');
        
        MortalityRandomGeneratorMale   = @(Stage, R)(uint16(floor(FmortMale(Stage,R))));
        MortalityRandomGeneratorFemale = @(Stage, R)(uint16(floor(FmortFemale(Stage,R))));
        
        %stage random number generator definition
        ageGroups = uint8([0 (19:5:85)+1]);
        ageGroup = sum(currentYear >= ageGroups); %which age group we are considering now
        
        femalesIndcs = cellfun(@(x)(strcmpi(x,'female')),stats.stageDistribution.sex);
        dataF = table2array(stats.stageDistribution(femalesIndcs,3:6));
        dataM = table2array(stats.stageDistribution(~femalesIndcs,3:6));
        dataF = dataF(ageGroup,:); dataFcum = [0 cumsum(dataF)];
        dataM = dataM(ageGroup,:); dataMcum = [0 cumsum(dataM)];
        StageProbability = [dataMcum; dataFcum]';
        
        meshStageCDF = unique(StageProbability);
        StageInvCDF = zeros(length(meshStageCDF),2);
        for i = 1:2
            StageInvCDF(:,i) = interp1(StageProbability(:,i)', 7:11, meshStageCDF);
        end
        
        %Gender == 2 is female, 1 is male
        Fstage = griddedInterpolant({1:2, meshStageCDF},StageInvCDF');
        StageRandomGenerator = @(Gender, R)(floor(Fstage(Gender,R)));
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  people die of natural causes     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diedNaturalCauses = rand(NAlive,1) < LifeTable(Gender,currentYear) & ~BlockFromDeath;%those will die now from other causes
    if any(diedNaturalCauses)
        %DeathYear(SubjectIDs(diedNaturalCauses)) = stepCounter;
        death(diedNaturalCauses, false);
    end
    
    %calculations for those that died before not from natural causes
    %(Included substitute)
    if ~isempty(GenderWouldBeAlive)
        diedNaturalCausesNotIncluded = rand(length(GenderWouldBeAlive),1) < LifeTable(GenderWouldBeAlive,currentYear);%those will die now from other causes
        if any(diedNaturalCausesNotIncluded)
            GenderWouldBeAlive(diedNaturalCausesNotIncluded) = [];
            SubjectsIDWouldBeAlive(diedNaturalCausesNotIncluded) = [];
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    people die of cancer           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diedOfCancer = Detected.MortTime < ni*numPeriods+1 & (stepCounter - Detected.CancerYear) >= Detected.MortTime;
    if any(diedOfCancer)
        deathIDs = unique(Detected.SubjectID(diedOfCancer));
        %DeathYear(deathIDs) = stepCounter;
        death(ismember(SubjectIDs, deathIDs));
    end
    
    R = rand(NAlive,2); %generating random numbers in advance
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a NEW POLYP appears               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    developNewPolyp = R(:,1) < PolypRate; %those will get new polyps
    numNewPolyps = sum(developNewPolyp);
    if numNewPolyps %if there is a new polyp to add; there is no limit on the number of polyps
        Polyp.SubjectID = [Polyp.SubjectID; SubjectIDs(developNewPolyp)];
        Polyp.Polyps = [Polyp.Polyps; ones(numNewPolyps,1)];
        Polyp.PolypYear = [Polyp.PolypYear; repmat(stepCounter,numNewPolyps,1)];
        location = LocationMatrix(randi(1000,numNewPolyps,1));
        Polyp.PolypLocation = [Polyp.PolypLocation; location];
        Polyp.OwnersGender = [Polyp.OwnersGender; Gender(developNewPolyp)];
        
        progressionRate = randi(500,numNewPolyps,1);
        Polyp.EarlyProgression = [Polyp.EarlyProgression; progressionRate];
        if flag.Correlation
            Polyp.AdvProgression = [Polyp.AdvProgression; progressionRate];
        else
            Polyp.AdvProgression = [Polyp.AdvProgression; randi(500,numNewPolyps,1)];
        end
        
        %calculating initial progression probability
        Polyp.ProgressionProb = [Polyp.ProgressionProb; LocationProgression(location)...
            .* GenderProgression(1+(Gender(developNewPolyp) - 1)*GenderProgressionRows) ... %need to use linear indexing
            .* RiskDistribution.EarlyRisk(progressionRate)];
        
        aux = StageVariables.FastCancer(1).*LocationProgression(location,6).*...
            GenderProgression(6+(Gender(developNewPolyp) - 1)*GenderProgressionRows);
        if strcmp(DwellSpeed,'Fast')
            aux = aux.*RiskDistribution.EarlyRisk(progressionRate);
        end
        Polyp.DirectProgressionProb = [Polyp.DirectProgressionProb; aux];
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % a NEW Cancer appears DIRECTLY     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    developNewCancer = R(:,2) < DirectCancerRate(Gender,currentYear) * DirectCancerSpeed; %those will get new cancer
    if any(developNewCancer)
        
        numNewCancers = sum(developNewCancer);
        Ca.Cancer = [Ca.Cancer; 7*ones(numNewCancers,1)];
        Ca.CancerYear = [Ca.CancerYear; repmat(stepCounter,numNewCancers,1)];
        location = LocationMatrix(randi(1000,numNewCancers,1),2);
        Ca.CancerLocation = [Ca.CancerLocation; location];
        Ca.DwellTime = [Ca.DwellTime; zeros(numNewCancers,1)];
        IDs = SubjectIDs(developNewCancer);
        Ca.SubjectID = [Ca.SubjectID; IDs];
        
        %SympStage = StageRandomGenerator(rand(numNewCancers,1));
        SympStage = StageRandomGenerator(double(Gender(developNewCancer)), rand(numNewCancers,1));
        SympTimeAdd = SojournRandomGenerator(SympStage,rand(numNewCancers,1));
        
        SympTimeAddDouble = double(SympTimeAdd);
        SympStageShifted = SympStage-6;
        Ca.TimeStage_I = [Ca.TimeStage_I;  stepCounter + uint16(SympTimeAddDouble.*StageDurationC(SympStageShifted,1))];
        Ca.TimeStage_II = [Ca.TimeStage_II;  stepCounter + uint16(SympTimeAddDouble.*StageDurationC(SympStageShifted,2))];
        Ca.TimeStage_III = [Ca.TimeStage_III;  stepCounter + uint16(SympTimeAddDouble.*StageDurationC(SympStageShifted,3))];
        
        
        Ca.SympTime = [Ca.SympTime; stepCounter+SympTimeAdd];
        Ca.SympStage = [Ca.SympStage; SympStage];
        Ca.OwnersGender = [Ca.OwnersGender; Gender(developNewCancer)];
        
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %      a polyp progresses           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(Polyp.Polyps) %if there is any polyp
        
        progressionProbTotal = AgeProgression(Polyp.Polyps,currentYear).*Polyp.ProgressionProb;
        directProgressionProbTotal = AgeProgression(6,currentYear).*Polyp.DirectProgressionProb;
        %generating random numbers in advance
        Rin = rand(length(progressionProbTotal),2);
        
        progressing = Rin(:,1) < progressionProbTotal;
        progressingDirectlyToCancer = ~progressing & Rin(:,2) < directProgressionProbTotal;
        
        if any(progressing | progressingDirectlyToCancer)
            
            Polyp.Polyps(progressing) = Polyp.Polyps(progressing) + 1; %increase stage
            
            indirectProgression = (progressing & Polyp.Polyps > 6);
            progressionToCancer = indirectProgression | progressingDirectlyToCancer;
            
            %update progression probability for those polyps that progressed
            probUpdate = progressing & ~progressionToCancer;
            Polyp.ProgressionProb(probUpdate) = Polyp.ProgressionProb(probUpdate).*...
                LocationProgressionRatios((Polyp.Polyps(probUpdate)-1)*LocationProgressionRows+Polyp.PolypLocation(probUpdate)).*...
                GenderProgressionRatios(Polyp.Polyps(probUpdate)+(Polyp.OwnersGender(probUpdate) - 1)*GenderProgressionRows);
            
            Polyp.DirectProgressionProb(probUpdate) = Polyp.DirectProgressionProb(probUpdate).*...
                StageVariables.FastCancerRatios(Polyp.Polyps(probUpdate));
            
            advancedRiskSwitch = probUpdate & Polyp.Polyps == 5;
            if any(advancedRiskSwitch)
                aux = (RiskDistribution.AdvancedRisk(Polyp.AdvProgression(advancedRiskSwitch))./RiskDistribution.EarlyRisk(Polyp.EarlyProgression(advancedRiskSwitch)));
                Polyp.ProgressionProb(advancedRiskSwitch) = Polyp.ProgressionProb(advancedRiskSwitch).*aux;
                if strcmp(DwellSpeed,'Fast')
                    Polyp.DirectProgressionProb(advancedRiskSwitch) = Polyp.DirectProgressionProb(advancedRiskSwitch).*aux;
                end
            end
            
            if any(progressionToCancer) %any progressing to cancer
                numNewCancers = sum(progressionToCancer);
                Ca.Cancer = [Ca.Cancer; 7*ones(numNewCancers,1)];
                Ca.CancerYear = [Ca.CancerYear; repmat(stepCounter,numNewCancers,1)];
                Ca.CancerLocation = [Ca.CancerLocation; Polyp.PolypLocation(progressionToCancer)];
                Ca.DwellTime = [Ca.DwellTime; stepCounter - Polyp.PolypYear(progressionToCancer)];
                IDs = Polyp.SubjectID(progressionToCancer);
                Ca.SubjectID = [Ca.SubjectID; IDs];
                
                %SympStage = StageRandomGenerator(rand(numNewCancers,1));
                SympStage = StageRandomGenerator(double(Polyp.OwnersGender(progressionToCancer)),rand(numNewCancers,1));
                SympTimeAdd = SojournRandomGenerator(SympStage,rand(numNewCancers,1));
                
                SympTimeAddDouble = double(SympTimeAdd);
                SympStageShifted = SympStage-6;
                Ca.TimeStage_I = [Ca.TimeStage_I;  stepCounter + uint16(SympTimeAddDouble.*StageDurationC(SympStageShifted,1))];
                Ca.TimeStage_II = [Ca.TimeStage_II;  stepCounter + uint16(SympTimeAddDouble.*StageDurationC(SympStageShifted,2))];
                Ca.TimeStage_III = [Ca.TimeStage_III;  stepCounter + uint16(SympTimeAddDouble.*StageDurationC(SympStageShifted,3))];
                
                Ca.SympTime = [Ca.SympTime; stepCounter+SympTimeAdd];
                Ca.SympStage = [Ca.SympStage; SympStage];
                Ca.OwnersGender = [Ca.OwnersGender; Polyp.OwnersGender(progressionToCancer)];
                
                %deleting those polyps that are progressing to cancer
                deletePolyps(progressionToCancer);
            end
            
        end
        
    end %end if there are any polyps to progress
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %   a polyp shrinks or disappears      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numPolyps = length(Polyp.Polyps);
    if numPolyps
        healing = rand(numPolyps,1) < StageVariables.Healing(Polyp.Polyps);
        if any(healing)
            
            Polyp.Polyps(healing) = Polyp.Polyps(healing)-1;
            dissapear = Polyp.Polyps == 0; %those healed completely
            
            %update progression probability for those polyps that shrinked
            probUpdate = healing & ~dissapear;
            Polyp.ProgressionProb(probUpdate) = Polyp.ProgressionProb(probUpdate)./...
                LocationProgressionRatios(Polyp.Polyps(probUpdate)*LocationProgressionRows+Polyp.PolypLocation(probUpdate))./...
                GenderProgressionRatios(Polyp.Polyps(probUpdate)+1+(Polyp.OwnersGender(probUpdate) - 1)*GenderProgressionRows);
            
            Polyp.DirectProgressionProb(probUpdate) = Polyp.DirectProgressionProb(probUpdate)./...
                StageVariables.FastCancerRatios(Polyp.Polyps(probUpdate)+1);
            
            advancedRiskSwitch = probUpdate & Polyp.Polyps == 4;
            if any(advancedRiskSwitch)
                aux = (RiskDistribution.AdvancedRisk(Polyp.AdvProgression(advancedRiskSwitch))./RiskDistribution.EarlyRisk(Polyp.EarlyProgression(advancedRiskSwitch)));
                Polyp.ProgressionProb(advancedRiskSwitch) = Polyp.ProgressionProb(advancedRiskSwitch)./aux;
                if strcmp(DwellSpeed,'Fast')
                    Polyp.DirectProgressionProb(advancedRiskSwitch) = Polyp.DirectProgressionProb(advancedRiskSwitch)./aux;
                end
            end
            
            %delete polyps that dissapear
            if any(dissapear)
                deletePolyps(dissapear);
            end
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % symptom development               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    developedSymptoms = stepCounter >= Ca.SympTime; %those will develop symptoms
    if any(developedSymptoms)
        Colonoscopy(unique(Ca.SubjectID(developedSymptoms)),'Symp');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cancer Progression                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Stage_I_progression = Ca.Cancer == 7 & stepCounter >= Ca.TimeStage_I;
    Stage_II_progression = Ca.Cancer == 8 & stepCounter >= Ca.TimeStage_II;
    Stage_III_progression = Ca.Cancer == 9 & stepCounter >= Ca.TimeStage_III;
    if any(Stage_I_progression)
        Ca.Cancer(Stage_I_progression) = 8;
    end
    if any(Stage_II_progression)
        Ca.Cancer(Stage_II_progression) = 9;
    end
    if any(Stage_III_progression)
        Ca.Cancer(Stage_III_progression) = 10;
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % polyp and cancer surveillance     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if mod(stepCounter,numPeriods) == 1 %if we start a new year
        
        SurveillanceFlag = false(1,NAlive);
        dcurrentYear = double(currentYear);
        
        if flag.Polyp_Surveillance
            
            LastPolyp = Last.Polyp(SubjectIDs); %select for those that are alive
            LastAdvPolyp = Last.AdvPolyp(SubjectIDs); %select for those that are alive
            LastColonoscopy = Last.Colonoscopy(SubjectIDs); %select for those that are alive
            
            % a polyp removed 5 years ago if no colonoscopy
            % performed inbetween.
            SurveillanceFlag = SurveillanceFlag | (dcurrentYear - LastPolyp == 5 & dcurrentYear - LastColonoscopy >= 5);
            % between 5 and 9 years after polyp removal after last colonoscopy
            SurveillanceFlag = SurveillanceFlag | ((dcurrentYear - LastPolyp > 5 & dcurrentYear - LastPolyp <= 9) & dcurrentYear-LastColonoscopy >= 5);
            % an advanced polyp 3 years ago
            SurveillanceFlag = SurveillanceFlag | (dcurrentYear - LastAdvPolyp == 3 & dcurrentYear-LastColonoscopy >= 3);
            % 5 years intervals if an advanced polyp had been diagnosed
            SurveillanceFlag = SurveillanceFlag | (LastAdvPolyp ~= -100 & dcurrentYear - LastAdvPolyp >= 5 & dcurrentYear-LastColonoscopy >= 5);
            if flag.AllPolypFollowUp
                SurveillanceFlag = SurveillanceFlag | (LastPolyp ~= -100 & dcurrentYear - LastAdvPolyp >= 5 & dcurrentYear-LastColonoscopy >= 5);
            end
        end
        
        if flag.Cancer_Surveillance
           
            LastCancer = Last.Cancer(SubjectIDs);
            LastColonoscopy = Last.Colonoscopy(SubjectIDs); %select for those that are alive
            
            hadCancer = LastCancer ~= -100;
            % 1 year after cancer diagnosis
            SurveillanceFlag = SurveillanceFlag | (hadCancer & (dcurrentYear - LastCancer == 1 & dcurrentYear - LastColonoscopy == 1));
            % 4 years after cancer diagnosis
            SurveillanceFlag = SurveillanceFlag | (hadCancer & (dcurrentYear - LastCancer == 4 & dcurrentYear - LastColonoscopy == 3));
             % 5 years intervals after this
            SurveillanceFlag = SurveillanceFlag | (hadCancer & (dcurrentYear - LastCancer >= 5 & dcurrentYear - LastColonoscopy >= 5));
            
        end
        
        if any(SurveillanceFlag)
             Colonoscopy(SubjectIDs(SurveillanceFlag),'Foll');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    summarizing polyps             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % before a polyp can be detected and removed for screening we
        % summarize the prevalence of all polyps and cancers
        if ~isempty(Polyp.Polyps)
            A = [Polyp.SubjectID Polyp.Polyps Polyp.OwnersGender];
            A = sortrows(A);
            aux = find(diff([A(:,1); NumIndividualsToSimulate+1]));
            %idxIDs = unique(A(:,1));
            
            %NumPolyps(idxIDs, currentYear) = diff([0; aux]);
            tmpHelp = diff([0; aux]);
            PolypsSumm.ActNumPolyps(1,currentYear) = sum(tmpHelp > 0);
            PolypsSumm.ActNumPolyps(2,currentYear) = sum(tmpHelp > 1);
            PolypsSumm.ActNumPolyps(3,currentYear) = sum(tmpHelp > 2);
            PolypsSumm.ActNumPolyps(4,currentYear) = sum(tmpHelp > 3);
            PolypsSumm.ActNumPolyps(5,currentYear) = sum(tmpHelp > 4);
            
            
            %MaxPolyps(idxIDs, currentYear) = A(aux,2);
            tmpHelp = A(aux,2);
            tmpHelp2 = A(aux,3) == 1; %take the gender
            
            %AllPolyps(unique(Polyp.Polyps), currentYear) = diff([0; find(diff([sort(Polyp.Polyps); 12]))]);
            
            
            PolypsSumm.NumPolyps(1,currentYear)  = sum(tmpHelp>0);
            PolypsSumm.NumPolyps(5,currentYear)  = sum(tmpHelp>4);
            PolypsSumm.FracPolyps(:,currentYear) = PolypsSumm.NumPolyps(:,currentYear)/length(SubjectIDs)*100;
            
            PolypsSumm.NumPolypsMale(1,currentYear)  = sum(tmpHelp(tmpHelp2)>0);
            PolypsSumm.NumPolypsMale(5,currentYear)  = sum(tmpHelp(tmpHelp2)>4);
            PolypsSumm.FracPolypsMale(:,currentYear) = PolypsSumm.NumPolypsMale(:,currentYear)/sum(Gender == 1)*100;
            
            PolypsSumm.NumPolypsFemale(1,currentYear)  = sum(tmpHelp(~tmpHelp2)>0);
            PolypsSumm.NumPolypsFemale(5,currentYear)  = sum(tmpHelp(~tmpHelp2)>4);
            PolypsSumm.FracPolypsFemale(:,currentYear) = PolypsSumm.NumPolypsFemale(:,currentYear)/sum(Gender == 2)*100;
            
            IncidenceCounter(4,currentYear) = length(Gender);
            IncidenceCounter(5,currentYear) = sum(Gender == 1);
            IncidenceCounter(6,currentYear) = sum(Gender == 2);
            
        end
        
        
        %fprintf('Calculating year %d\n', currentYear);
    end
    
end

%% modyfing output arguments for compability with old CMOST
currentYear = double(currentYear);

%TumorRecord.Time = double(TumorRecord.Time)/double(numPeriods) + 0.75;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS DEFINITIONS            %
% it is crucial to have them nested, i.e. defined before the last end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function death(whichSubjects, keep)
        if nargin < 2
            keep = true; %default value is to keep for further calculations
        end
        %removing all associated polyps
        IDsToDelete = SubjectIDs(whichSubjects);
        indx = ismember(Polyp.SubjectID, IDsToDelete);
        deletePolyps(indx);
        
        %removing all associated cancers
        indx = ismember(Ca.SubjectID, IDsToDelete);
        deleteCancers(indx);
        
        %removing all associated detected cancers
        indx = ismember(Detected.SubjectID, IDsToDelete);
        deleteDetectedCancers(indx);
        
        if keep %if too keep individuals for further, Included substitute
            SubjectsIDWouldBeAlive = [SubjectsIDWouldBeAlive; SubjectIDs(whichSubjects)];
            GenderWouldBeAlive     = [GenderWouldBeAlive; Gender(whichSubjects)];
        end
        
        %making variable consitent with CMOST old codes
        %for ii = 1:length(IDsToDelete)
        %    NumCancer(IDsToDelete(ii), currentYear:end) = NumCancer(IDsToDelete(ii), currentYear);
        %    MaxCancer(IDsToDelete(ii), currentYear:end) = MaxCancer(IDsToDelete(ii), currentYear);
        %end
        
        %removing individuals from simulation when they die
        whichSubjects = ~whichSubjects;
        
        SubjectIDs = SubjectIDs(whichSubjects);
        Gender = Gender(whichSubjects);
        NAlive = sum(whichSubjects);
        IndividualRisk = IndividualRisk(whichSubjects);
        PolypRate = PolypRate(whichSubjects);
        BlockFromDeath = BlockFromDeath(whichSubjects);
        
        ScreeningPreference = ScreeningPreference(whichSubjects);
    end

    function deletePolyps(indx)
        
        if ~islogical(indx)
            indxN = true(size(Polyp.SubjectID));
            indxN(indx) = false;
            indx = indxN;
        else
            indx = ~indx;
        end
        
        Polyp.SubjectID = Polyp.SubjectID(indx);
        Polyp.Polyps = Polyp.Polyps(indx);
        Polyp.PolypYear = Polyp.PolypYear(indx);
        Polyp.PolypLocation = Polyp.PolypLocation(indx);
        Polyp.EarlyProgression = Polyp.EarlyProgression(indx);
        Polyp.AdvProgression = Polyp.AdvProgression(indx);
        Polyp.ProgressionProb = Polyp.ProgressionProb(indx);
        Polyp.DirectProgressionProb = Polyp.DirectProgressionProb(indx);
        Polyp.OwnersGender = Polyp.OwnersGender(indx);
        %         end
    end

    function deleteCancers(indx)
        
        if ~islogical(indx)
            indxN = true(size(Ca.Cancer));
            indxN(indx) = false;
            indx = indxN;
        else
            indx = ~indx;
        end
        
        Ca.Cancer               = Ca.Cancer(indx);
        Ca.CancerYear           = Ca.CancerYear(indx);
        Ca.CancerLocation       = Ca.CancerLocation(indx);
        Ca.DwellTime            = Ca.DwellTime(indx);
        Ca.SympTime             = Ca.SympTime(indx);
        Ca.SympStage            = Ca.SympStage(indx);
        Ca.TimeStage_I          = Ca.TimeStage_I(indx);
        Ca.TimeStage_II         = Ca.TimeStage_II(indx);
        Ca.TimeStage_III        = Ca.TimeStage_III(indx);
        Ca.SubjectID            = Ca.SubjectID(indx);
        Ca.OwnersGender         = Ca.OwnersGender(indx);
    end

    function deleteDetectedCancers(indx)
        if ~islogical(indx)
            indxN = true(size(Detected.Cancer));
            indxN(indx) = false;
            indx = indxN;
        else
            indx = ~indx;
        end
        
        Detected.Cancer               = Detected.Cancer(indx);
        Detected.CancerYear           = Detected.CancerYear(indx);
        Detected.CancerLocation       = Detected.CancerLocation(indx);
        Detected.MortTime             = Detected.MortTime(indx);
        Detected.SubjectID            = Detected.SubjectID(indx);
    end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%         COLONOSCOPY                               %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function Colonoscopy(SelectedSubjectID, Modus)
        
        %check whether it is PBP mode
        PBPudpate = any(strcmpi(Modus,{'PBPx','PBPm'}));
        
        Ncolo = length(SelectedSubjectID); %number of colonoscopies to perform
        
        [ScreenedPolyps, LocP] = ismember(Polyp.SubjectID, SelectedSubjectID);
        [ScreenedCancers, LocC] = ismember(Ca.SubjectID, SelectedSubjectID);
        LocP(LocP==0) = 1; %taking care of zero indexing that would throw an error
        LocC(LocC==0) = 1; %taking care of zero indexing that would throw an error
        
        %select only those polyps and cancers that are reached by
        %colonoscopy (cecum = 1, rectum = 13))
        CurrentReach = ColoReachRandomGenerator(rand(Ncolo,1));
        ScreenedPolyps = ScreenedPolyps & (Polyp.PolypLocation >= CurrentReach(LocP));
        ScreenedCancers = ScreenedCancers & (Ca.CancerLocation >= CurrentReach(LocC));
        
        nScreenedPolyps = sum(ScreenedPolyps);
        nScreenedCancers = sum(ScreenedCancers);
        
        Fpolyps = find(ScreenedPolyps);
        detectedPolyps = rand(nScreenedPolyps,1) < StageVariables.Colo_Detection(Polyp.Polyps(ScreenedPolyps)).*...
            Location.ColoDetection(Polyp.PolypLocation(ScreenedPolyps));
        
        
        counterAdvanced = uint8(zeros(Ncolo,1)); %number of detected advanced polyps per patient
        counterEarly    = uint8(zeros(Ncolo,1)); %number of detected early polyps per patient
        
        if any(detectedPolyps)
            idx = Fpolyps(detectedPolyps);
            X = LocP(idx);
            advanced = Polyp.Polyps(idx)>4;
            U = unique(X);
            for ii = 1:length(U)
                counterAdvanced(U(ii)) = sum(X == U(ii) & advanced);
                counterEarly(U(ii)) = sum(X == U(ii) & ~advanced);
            end
            
            Last.Polyp(SelectedSubjectID(counterEarly > 0)) = currentYear;
            Last.AdvPolyp(SelectedSubjectID(counterAdvanced > 0 | counterEarly > 2)) = currentYear;% 3 polyps counts as an advanced polyp
            
            % in case of detecion we remove the polyp
            deletePolyps(idx); %deleting detected polyps
        end
        
      
        
        [~, m] = ismember(Modus,{'Scre','Symp','Foll','Base'});
        m = uint8(m);
        if PBPudpate
            m = uint8(1);
        end
        
        Fcancers = find(ScreenedCancers);
        detectedCancers = rand(nScreenedCancers,1) < StageVariables.Colo_Detection(Ca.Cancer(ScreenedCancers));
        
        counterCancer    = uint8(zeros(Ncolo,1)); %number of detected cancers per patient
        
        if any(detectedCancers)
            % in case of detecion we mark cancer as detected and remove it
            idx = Fcancers(detectedCancers);
            
            X = LocC(idx);
            U = unique(X);
            for ii = 1:length(U)
                counterCancer(U(ii)) = sum(X == U(ii));
            end
            
            Detected.Cancer         = [Detected.Cancer; Ca.Cancer(idx)];
            Detected.CancerYear     = [Detected.CancerYear; repmat(stepCounter, length(idx),1)];
            Detected.CancerLocation = [Detected.CancerLocation; Ca.CancerLocation(idx)];
            %MortTime = MortalityRandomGenerator(double(Ca.Cancer(idx)),rand(length(idx),1));%+4;
            MortTime = zeros(length(idx),1);
            malesLoc = Ca.OwnersGender(idx) == 1;
            MortTime(malesLoc) = MortalityRandomGeneratorMale(double(Ca.Cancer(idx(malesLoc))),rand(sum(malesLoc),1));%+4;
            MortTime(~malesLoc) = MortalityRandomGeneratorFemale(double(Ca.Cancer(idx(~malesLoc))),rand(sum(~malesLoc),1));%+4;
            
            %we need to block the patient from dying before the MortTime
            %first I need to locate the patients in the whole cohort
            [im, LocC] = ismember(unique(Ca.SubjectID(idx(MortTime < ni*numPeriods+1))), SubjectIDs);
            if ~all(im)
                disp('Grave error #ID colonoscopyBlockPatient');
            else
                %finding which patients to block
                %disp('Blocking')
                
                BlockFromDeath(LocC) = true;
            end
            
            Detected.MortTime       = [Detected.MortTime; MortTime];
            Detected.SubjectID      = [Detected.SubjectID; Ca.SubjectID(idx)];
            
            % we need keep track of key parameters
            % Stage Location SojournTime Sex DetectionMode
            X = Ca.SubjectID(idx);%IDs of the patients with detected cancers
            U = unique(X);
            N = sum(counterCancer > 0);
            idxFin = zeros(N,1);
            for ii = 1:length(U) %for each patient find the one with the most advanced stage
                idxF = idx(X == U(ii));
                [~, whichMax] = max(Ca.Cancer(idxF));
                idxFin(ii) = idxF(whichMax(end)); %this cancer will be saved
            end
            
            %here we update incidence
            IncidenceCounter(1,currentYear) = IncidenceCounter(1,currentYear) + N; %overall
            IncidenceCounter(2,currentYear) = IncidenceCounter(2,currentYear) + sum(Ca.OwnersGender(idxFin)==1); %males
            IncidenceCounter(3,currentYear) = IncidenceCounter(3,currentYear) + sum(Ca.OwnersGender(idxFin)==2); %females
            
            Last.Cancer(SelectedSubjectID(counterCancer > 0)) = currentYear;
            
            deleteCancers(idx); %deleting deetcted cancers from the array
        end

        
        Last.Colonoscopy(SelectedSubjectID) = currentYear;
        
        factor = 1.5*ones(Ncolo,1);
        %noTumorAndNoPolyps = counterCancer == 0 & counterAdvanced == 0 & counterEarly == 0;
        
        %TumorAndNoPolyps = counterCancer > 0  & counterAdvanced == 0 & counterEarly == 0;
        
        %%%% Complications
        Rand = rand(Ncolo,6); %generating random numbers in advance
        
        perforation  = Rand(:,1) < risc.Colonoscopy_RiscPerforation.*factor;
        if any(perforation)% a perforation happend
            
            deathPerforation        = perforation & Rand(:,2) < risc.DeathPerforation;
            if any(deathPerforation) % patient died during colonoscopy from a perforation
                %DeathYear(SelectedSubjectID(deathPerforation)) = stepCounter;
                
                death(ismember(SubjectIDs, SelectedSubjectID(deathPerforation)));
            end
        end
        serosaBurn = ~perforation & Rand(:,3) < risc.Colonoscopy_RiscSerosaBurn.*factor;
        
        bleeding = ~serosaBurn & ~perforation & Rand(:,4) < risc.Colonoscopy_RiscBleeding.*factor;
        
        bleedingTransfusion     = ~bleeding & ~serosaBurn & ~perforation  &  Rand(:,5) < risc.Colonoscopy_RiscBleedingTransfusion.*factor;
        if any(bleedingTransfusion)% bleeding recquiring transfusion
            
            deathBleedingTransfusion = bleedingTransfusion & Rand(:,6) < risc.DeathBleedingTransfusion;
            if any(deathBleedingTransfusion) % patient died during colonoscopy from a bleeding complication
                %DeathYear(SelectedSubjectID(deathBleedingTransfusion)) = stepCounter;
                
                death(ismember(SubjectIDs, SelectedSubjectID(deathBleedingTransfusion)));
            end
            
        end
        
        
    end


end


