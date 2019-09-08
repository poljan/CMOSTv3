
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     CMOST: Colon Modeling with Open Source Tool
%     created by Meher Prakash and Benjamin Misselwitz 2012 - 2016
%
%     This program is part of free software package CMOST for colo-rectal  
%     cancer simulations: You can redistribute it and/or modify 
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%       
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EvaluatePBP_23_04_2018

ExcelFileName              = 'PBP_CMOST13_23042018';
ExcelSavePath              = '/Users/misselwb/Documents/CRC Screening/_Bowel Prep_2/Data';
ExcelFileName              = fullfile(ExcelSavePath, ExcelFileName);

AnalysisPipeline = '/Users/misselwb/Documents/CRC Screening/_Bowel Prep_2/Data/PBP';
cd (AnalysisPipeline)
MasterDirectory=dir;
for f=1:length(MasterDirectory)
    AllFiles{f} = MasterDirectory(f).name;
end
z=0;
Error = zeros(length(MasterDirectory), 1);
for x=1:length(MasterDirectory)
    if ~or(or(strcmp(MasterDirectory(x,1).name, '.'), strcmp(MasterDirectory(x,1).name, '..')) ,...
            or(strcmp(MasterDirectory(x,1).name, '.DS_S_Results'), strcmp(MasterDirectory(x,1).name, '.DS_S_Store')))
        if and(~MasterDirectory(x,1).isdir, isempty(regexp(MasterDirectory(x,1).name, 'Results.mat$', 'once')))
            z=z+1;
            l = length(MasterDirectory(x,1).name)-4;
            FileName{z}  = [MasterDirectory(x,1).name(1:l) '_Results'];
            Directory{z} = pwd;
            if ~isempty(find(ismember(AllFiles, [FileName{z} '.mat'])==1))
                Error(z) = 0;
            else
                Error(z) = 1;
            end
        end
    end
end
if max(Error) > 0
    display('Error, for the following files no output file was detected')
    for f=1:length(Error)
        if Error(f)
            display(Directory{f})
            display(FileName{f})
        end
    end
    return
end
BatchCounter = z;

HeadLineColumns=...
    {'B2' 'C2'  'D2'  'E2'  'F2'  'G2'  'H2'  'I2'  'J2'  'K2'  'L2'  'M2'  'N2'  'O2'  'P2'  'Q2' 'R2' 'S2' 'T2'...
    'U2'  'V2'  'W2'  'X2'  'Y2'  'Z2'  'AA2' 'AB2' 'AC2' 'AD2' 'AE2' 'AF2' 'AG2' 'AH2' 'AI2' 'AJ2'...
    'AK2' 'AL2' 'AM2' 'AN2' 'AO2' 'AP2' 'AQ2' 'AR2' 'AS2' 'AT2' 'AU2' 'AV2' 'AW2' 'AX2' 'AY2' 'AZ2'...
    'BA2' 'BB2' 'BC2' 'BD2' 'BE2' 'BF2' 'BG2' 'BH2' 'BI2' 'BJ2' 'BK2' 'BL2' 'BM2' 'BN2'...
    'Bo2' 'Bp2' 'Bq2' 'Br2' 'Bs2' 'Bt2' 'Bu2' 'Bv2' 'Bw2' 'Bx2' 'By2' 'Bz2' 'ca2' 'cb2'...
    'CC2' 'CD2' 'CE2' 'CF2' 'CG2' 'CH2' 'CI2' 'CJ2' 'CK2' 'CL2' 'CM2' 'CN2'...
    'CO2' 'CP2' 'CQ2' 'CR2' 'CS2' 'CT2' 'CU2' 'CV2' 'CW2' 'CX2' 'CY2' 'CZ2' 'DA2' 'DB2'...
    'DC2' 'DD2' 'DE2' 'DF2' 'DG2' 'DH2' 'DI2' 'DJ2' 'DK2' 'DL2' 'DM2' 'DN2'};

clear Headlines
Headlines{1,1} =  'no_intervention';
Headlines{1,2} =  'Screening_Kolo';
Headlines{1,3} =  '_SC01_Aronchick_01_avg_';
Headlines{1,4} =  '_SC01_Aronchick_01_low_';
Headlines{1,5} =  '_SC01_Aronchick_01_high_';
Headlines{1,6} =  '_SC01_Aronchick_02_avg_';
Headlines{1,7} =  '_SC01_Aronchick_02_low_';
Headlines{1,8} =  '_SC01_Aronchick_02_high_';
Headlines{1,9} =  '_SC01_Aronchick_03_avg_';
Headlines{1,10} = '_SC01_Aronchick_03_low_';
Headlines{1,11} = '_SC01_Aronchick_03_high_';
Headlines{1,12} = '_SC01_Aronchick_04_avg_';
Headlines{1,13} = '_SC01_Aronchick_04_low_';
Headlines{1,14} = '_SC01_Aronchick_04_high_';

Headlines{1,15} =  '_SC02_PBP_66_WWxx_Aronchick_01_avg_'; % insufficient
Headlines{1,16} =  '_SC02_PBP_66_WWxx_Aronchick_02_avg_'; % poor
Headlines{1,17} =  '_SC02_PBP_66_WWxx_Aronchick_03_avg_'; % fair
Headlines{1,18} =  '_SC02_PBP_66_WWxx_Aronchick_04_avg_'; % good
Headlines{1,19} =  '_SC02_PBP_66_WWxx_Aronchick_05_avg_'; % excellent

Headlines{1,20} = '_SC03_PBP_66_WW00_Aronchick_01_avg_';
Headlines{1,21} = '_SC03_PBP_66_WW01_Aronchick_01_avg_';
Headlines{1,22} = '_SC03_PBP_66_WW02_Aronchick_01_avg_';
Headlines{1,23} = '_SC03_PBP_66_WW03_Aronchick_01_avg_';
Headlines{1,24} = '_SC03_PBP_66_WW04_Aronchick_01_avg_';
Headlines{1,25} = '_SC03_PBP_66_WW05_Aronchick_01_avg_';
Headlines{1,26} = '_SC03_PBP_66_WW06_Aronchick_01_avg_';
Headlines{1,27} = '_SC03_PBP_66_WW07_Aronchick_01_avg_';
Headlines{1,28} = '_SC03_PBP_66_WW08_Aronchick_01_avg_';
Headlines{1,29} = '_SC03_PBP_66_WW09_Aronchick_01_avg_';
Headlines{1,30} = '_SC03_PBP_66_WW10_Aronchick_01_avg_';

Headlines{1,31} = '_SC03_PBP_66_WW00_Aronchick_02_avg_';
Headlines{1,32} = '_SC03_PBP_66_WW01_Aronchick_02_avg_';
Headlines{1,33} = '_SC03_PBP_66_WW02_Aronchick_02_avg_';
Headlines{1,34} = '_SC03_PBP_66_WW03_Aronchick_02_avg_';
Headlines{1,35} = '_SC03_PBP_66_WW04_Aronchick_02_avg_';
Headlines{1,36} = '_SC03_PBP_66_WW05_Aronchick_02_avg_';
Headlines{1,37} = '_SC03_PBP_66_WW06_Aronchick_02_avg_';
Headlines{1,38} = '_SC03_PBP_66_WW07_Aronchick_02_avg_';
Headlines{1,39} = '_SC03_PBP_66_WW08_Aronchick_02_avg_';
Headlines{1,40} = '_SC03_PBP_66_WW09_Aronchick_02_avg_';
Headlines{1,41} = '_SC03_PBP_66_WW10_Aronchick_02_avg_';

Headlines{1,42} = '_SC03_PBP_66_WW00_Aronchick_03_avg_';
Headlines{1,43} = '_SC03_PBP_66_WW01_Aronchick_03_avg_';
Headlines{1,44} = '_SC03_PBP_66_WW02_Aronchick_03_avg_';
Headlines{1,45} = '_SC03_PBP_66_WW03_Aronchick_03_avg_';
Headlines{1,46} = '_SC03_PBP_66_WW04_Aronchick_03_avg_';
Headlines{1,47} = '_SC03_PBP_66_WW05_Aronchick_03_avg_';
Headlines{1,48} = '_SC03_PBP_66_WW06_Aronchick_03_avg_';
Headlines{1,49} = '_SC03_PBP_66_WW07_Aronchick_03_avg_';
Headlines{1,50} = '_SC03_PBP_66_WW08_Aronchick_03_avg_';
Headlines{1,51} = '_SC03_PBP_66_WW09_Aronchick_03_avg_';
Headlines{1,52} = '_SC03_PBP_66_WW10_Aronchick_03_avg_';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Scenario 1                         %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Screening %%%
% discounting after year 50
DisCountMask = zeros(101,1);
DisCountMask(1:51)=1;
for ff=52:101
    DisCountMask(ff) = DisCountMask(ff-1)* 0.97;
end

BatchCounter = z;
NoInterventionCounter =0;
StandardColoCounter   =0;
ScreenCounter = zeros(1, 10);

for z=1:BatchCounter
    pos1  = regexp(FileName{z}, '_SC01_Baseline_');
    pos2  = regexp(FileName{z}, '_SC01_Standard_Screening_');
    
    flag = false;
    for f=3:length(Headlines)
        posX = regexp(FileName{z}, Headlines{z}, 'once');
        if ~isempty(posX)
            flag = true;
            pos = f;
        end
    end

    if ~isempty(pos1) || ~isempty(pos2) || flag

        clear Results
        load(fullfile(Directory{z}, FileName{z}));
        
        Results.TotalCosts   = Results.Treatment + Results.Screening + Results.FollowUp + Results.Other;
        Results.DiscCosts    = Results.TotalCosts .* DisCountMask(1:100);
        try
            Results.SumDiscCosts = sum(Results.DiscCosts(51:100));
            
            Results.YearsLost    = Results.YearsLostCa + Results.YearsLostColo;
            Results.DiscYears    = Results.YearsLost .* DisCountMask(1:length(Results.YearsLost))';
            Results.SumDiscYears = sum(Results.DiscYears(51:100));
        catch exception
            rethrow(exception)
        end
        if ~isempty(pos1) % no intervention
            num = 1; ScreenCounter(num) = ScreenCounter(num)+1;
            AllResults{num, ScreenCounter(num)} = Results.Variable;
            NoInterventionCounter = NoInterventionCounter +1;
            control.NoIntervention.CaDeath(NoInterventionCounter)      = Results.Variable{13};
            control.NoIntervention.LYlost(NoInterventionCounter)       = Results.Variable{14};
            control.NoIntervention.AllCa(NoInterventionCounter)        = Results.Variable{37};
            control.NoIntervention.YLostColo(NoInterventionCounter)    = Results.Variable{16};
            control.NoIntervention.TotalCosts(NoInterventionCounter)   = Results.Variable{17};
            control.NoIntervention.DeathColo(NoInterventionCounter)    = Results.Variable{15};
            control.NoIntervention.SumDiscYears(NoInterventionCounter) = Results.SumDiscYears;
            control.NoIntervention.SumDiscCosts(NoInterventionCounter) = Results.SumDiscCosts;
            NumPat(num, ScreenCounter(num))  = Results.NumberPatients(51);
        elseif ~isempty(pos2) % screening colo
            num = 2; ScreenCounter(num) = ScreenCounter(num)+1;
            AllResults{num, ScreenCounter(num)} = Results.Variable;
            StandardColoCounter = StandardColoCounter +1;
            control.StandardColo.CaDeath(StandardColoCounter)      = Results.Variable{13};
            control.StandardColo.LYlost(StandardColoCounter)       = Results.Variable{14};
            control.StandardColo.AllCa(StandardColoCounter)        = Results.Variable{37};
            control.StandardColo.YLostColo(StandardColoCounter)    = Results.Variable{16};
            control.StandardColo.TotalCosts(StandardColoCounter)   = Results.Variable{17};
            control.StandardColo.DeathColo(StandardColoCounter)    = Results.Variable{15};
            control.StandardColo.SumDiscYears(StandardColoCounter) = Results.SumDiscYears;
            control.StandardColo.SumDiscCosts(StandardColoCounter) = Results.SumDiscCosts;
            NumPat(num, ScreenCounter(num))  = Results.NumberPatients(51);
        else 
            ScreenCounter(pos) = ScreenCounter(pos)+1;
            AllResults{pos, ScreenCounter(pos)} = Results.Variable;
            NumPat(pos, ScreenCounter(pos))     = Results.NumberPatients(51);
        end
    end
end

for f=1:length(Headlines)
    AllResultsVector{f} = cell(67,1);
end

for x1=1:length(Headlines) % we re-assemble everything 
    for x2=1:ScreenCounter(x1)
        for x3=1:length(AllResults{x1, x2})
            if isnumeric(AllResults{x1, x2}{x3})
                AllResultsVector{x1}{x3}(x2) = AllResults{x1, x2}{x3};
            end
        end
    end
end

NumPatAll = mean(mean(NumPat));
SummaryLegend = Results.Var_Legend;

% incidene reduction
SummaryLegend{60} = 'incidence reduction';
for f=1:length(Headlines)
    AllResultsVector{f}{60} = (AllResultsVector{1}{37}-AllResultsVector{f}{37}) ./AllResultsVector{1}{37}*100;
end

% mortality reduction
SummaryLegend{61} = 'mortality reduction';
for f=1:length(Headlines)
    AllResultsVector{f}{61} = (AllResultsVector{1}{13}-AllResultsVector{f}{13}) ./AllResultsVector{1}{13}*100;
end

% life years gained undiscounted
SummaryLegend{62} = 'life years gained';
for f=1:length(Headlines)
    AllResultsVector{f}{62} = ((AllResultsVector{1}{14}+AllResultsVector{1}{16}) - (AllResultsVector{f}{14}+AllResultsVector{f}{16})) /100;
end

% screening colonoscopies
SummaryLegend{63} = 'screening colonoscopies';
for f=1:length(Headlines)
    AllResultsVector{f}{63} = AllResultsVector{f}{5}/100;
end

% surveillance colonoscopies
SummaryLegend{64} = 'surveillance colonoscopies';
for f=1:length(Headlines)
    AllResultsVector{f}{64} = AllResultsVector{f}{7}/100;
end

% all colonoscopies
SummaryLegend{65} = 'all colonoscopies';
for f=1:length(Headlines)
    AllResultsVector{f}{65} = (AllResultsVector{f}{5}+AllResultsVector{f}{6}+AllResultsVector{f}{7})/100;
end

% CRC cases prevented
SummaryLegend{66} = 'crc cases prevented';
for f=1:length(Headlines)
    AllResultsVector{f}{66} = (AllResultsVector{1}{37}-AllResultsVector{f}{37})/100;
end

% CRC mortality cases prevented
SummaryLegend{67} = 'crc mortality cases prevented';
for f=1:length(Headlines)
    AllResultsVector{f}{67} = (AllResultsVector{1}{13}-AllResultsVector{f}{13})/100;
end

% colonoscopies per case prevented
SummaryLegend{68} = 'colonoscopies per case prevented';
for f=1:length(Headlines)
    AllResultsVector{f}{68} = AllResultsVector{f}{65}./AllResultsVector{f}{66};
end

% colonoscopies per life year gained
SummaryLegend{69} = 'colonoscopies per LY gained';
try
for f=1:length(Headlines)
    AllResultsVector{f}{69} = AllResultsVector{f}{65}./AllResultsVector{f}{62};
end
catch exception
    rethrow(exception)
end

% we calculate 
for f=1:length(Headlines) % we calculate an average
    for x3=1:length(AllResultsVector{f, 1})
        AllResultsMean{f}{x3} = mean(AllResultsVector{f}{x3});
        AllResultsStd{f}{x3}  = std(AllResultsVector{f}{x3});
    end
end


xlswrite (ExcelFileName, Headlines, 'PBP_mean', 'B1')
xlswrite (ExcelFileName, Headlines, 'PBP_std', 'B1')
xlswrite (ExcelFileName, SummaryLegend, 'PBP_mean', 'A2')
xlswrite (ExcelFileName, SummaryLegend, 'PBP_std', 'A2')
for f=1:length(Headlines)
    xlswrite (ExcelFileName, AllResultsMean{f}, 'PBP_mean', HeadLineColumns{f})
    xlswrite (ExcelFileName, AllResultsStd{f},  'PBP_std',  HeadLineColumns{f})
end

%%%%%%%%%%%%%%%%%%%%%%
%%%  Scenario 2    %%%
%%%%%%%%%%%%%%%%%%%%%%

ReplicaNumber = 0;

for z=1:BatchCounter 
    flag = false;
    for f=15:length(Headlines)
        posX = regexp(FileName{z}, Headlines{z}, 'once'); 
        if ~isempty(posX)
            flag = true;
            pos  = posX;
        end
    end
        
    if flag 
        FileTemp = FileName{z};
        posX = regexp(FileTemp, '_'); posX = posX(end);
        NumTmp1= str2num(FileTemp((posX+1)));
        NumTmp2= str2num(FileTemp((posX+2)));
        NumTmp3= str2num(FileTemp((posX+3)));
        if isempty(NumTmp2)
            ReplicaVar = NumTmp1;
        elseif isempty(NumTmp3)
            ReplicaVar = NumTmp1*10+NumTmp2;
        else
            ReplicaVar = NumTmp1*100+NumTmp1*10+NumTmp2;
        end
        ReplicaNumber = max(ReplicaNumber, ReplicaVar);
        
        clear Results 
        Results = load(fullfile(Directory{z}, FileTemp));
        
        % DiagnosedCancer(y, z); the number indicates the maximum stage
        % timor diagnosed at that year
        
        % all individuals
        DiagnosedCancer = Results.DiagnosedCancer;
        AliveArea = zeros(n, 100);
        
        % we loop over all individuals and see whether they were alive
        for f=1:n
            if Results.DeathYear(f) <= 66 
                DiagnosedCancer(:, f) = 0;
            else
                AliveArea(f, 1 : round(Results.DeathYear(f))) = 1;
            end
        end
        IntervallCancer.All{pos, ReplicaVar} = sum(DiagnosedCancer, 2);
        Alive.All{pos, ReplicaVar} = sum(AliveArea, 1);
        
        % NO POLYP
        if pos >14
            DiagnosedCancer = Results.DiagnosedCancer;
            AliveArea = zeros(n, 100);
            
            % we loop over all individuals and see whether they were alive and
            % have no polyp
            for f=1:n
                if Results.DeathYear(f) > 66 && ~isequal(Results.PBC_Doc.Early, 0)...
                        && ~isequal(Results.PBC_Doc.Advanced, 0) && ~isequal(Results.PBC_Doc.Cancer, 0)
                    AliveArea(f, 1 : round(Results.DeathYear(f))) = 1;
                else
                    DiagnosedCancer(:, f) = 0;
                end
            end
            IntervallCancer.NoPolyp{pos, ReplicaVar} = sum(DiagnosedCancer, 2);
            Alive.NoPolyp{pos, ReplicaVar} = sum(AliveArea, 1);
        end
        
        % ONE EARLY POLYP
        if pos >14
            DiagnosedCancer = Results.DiagnosedCancer;
            AliveArea = zeros(n, 100);
            
            % we loop over all individuals and see whether they were alive and
            % have no polyp
            for f=1:n
                if Results.DeathYear(f) > 66 && ~isequal(Results.PBC_Doc.Early, 1)...
                        && ~isequal(Results.PBC_Doc.Advanced, 0) && ~isequal(Results.PBC_Doc.Cancer, 0)
                    AliveArea(f, 1 : round(Results.DeathYear(f))) = 1;
                else
                    DiagnosedCancer(:, f) = 0;
                end
            end
            IntervallCancer.OnePolyp{pos, ReplicaVar} = sum(DiagnosedCancer, 2);
            Alive.OnePolyp{pos, ReplicaVar} = sum(AliveArea, 1);
        end
        
        % TWO EARLY POLYPs
        if pos >14
            DiagnosedCancer = Results.DiagnosedCancer;
            AliveArea = zeros(n, 100);
            
            % we loop over all individuals and see whether they were alive and
            % have no polyp
            for f=1:n
                if Results.DeathYear(f) > 66 && ~isequal(Results.PBC_Doc.Early, 2)...
                        && ~isequal(Results.PBC_Doc.Advanced, 0) && ~isequal(Results.PBC_Doc.Cancer, 0)
                    AliveArea(f, 1 : round(Results.DeathYear(f))) = 1;
                else
                    DiagnosedCancer(:, f) = 0;
                end
            end
            IntervallCancer.TwoPolyp{pos, ReplicaVar} = sum(DiagnosedCancer, 2);
            Alive.TwoPolyp{pos, ReplicaVar} = sum(AliveArea, 1);
        end
        
        % ONE ADVANCED POLYP
        if pos >14
            DiagnosedCancer = Results.DiagnosedCancer;
            AliveArea = zeros(n, 100);
            
            % we loop over all individuals and see whether they were alive and
            % have no polyp
            for f=1:n
                if Results.DeathYear(f) > 66 ...
                        && ~isequal(Results.PBC_Doc.Advanced, 1) && ~isequal(Results.PBC_Doc.Cancer, 0)
                    AliveArea(f, 1 : round(Results.DeathYear(f))) = 1;
                else
                    DiagnosedCancer(:, f) = 0;
                end
            end
            IntervallCancer.OneAdvPolyp{pos, ReplicaVar} = sum(DiagnosedCancer, 2);
            Alive.OneAdvPolyp{pos, ReplicaVar} = sum(AliveArea, 1);
        end
    end
end
        
for x1=15:length(Headlines) % we re-assemble everything
    for x2=1:ReplicaNumber
        if isequal(x2, 1)
            R_IntervallCancer.All{x1}         = IntervallCancer.All{x1, x2};
            R_Alive.All{x1}                   = Alive.All{x1, x2};
            R_IntervallCancer.NoPolyp{x1}     = IntervallCancer.NoPolyp{x1, x2};
            R_Alive.NoPolyp{x1}               = Alive.NoPolyp{x1, x2};
            R_IntervallCancer.OnePolyp{x1}    = IntervallCancer.OnePolyp{x1, x2};
            R_Alive.OnePolyp{x1}              = Alive.OnePolyp{x1, x2};
            R_IntervallCancer.TwoPolyp{x1}    = IntervallCancer.TwoPolyp{x1, x2};
            R_Alive.TwoPolyp{x1}              = Alive.TwoPolyp{x1, x2};
            R_IntervallCancer.OneAdvPolyp{x1} = IntervallCancer.OneAdvPolyp{x1, x2};
            R_Alive.OneAdvPolyp{x1}           = Alive.OneAdvPolyp{x1, x2};
        else
            R_IntervallCancer.All{x1}         = cat(2, R_IntervallCancer.All{x1}, IntervallCancer.All{x1, x2});
            R_Alive.All{x1}                   = cat(2, R_Alive.All{x1}, Alive.All{x1, x2});
            R_IntervallCancer.NoPolyp{x1}     = cat(2, R_IntervallCancer.NoPolyp{x1}, IntervallCancer.NoPolyp{x1, x2});
            R_Alive.NoPolyp{x1}               = cat(2, R_Alive.NoPolyp{x1}, Alive.NoPolyp{x1, x2});
            R_IntervallCancer.OnePolyp{x1}    = cat(2, R_IntervallCancer.OnePolyp{x1}, IntervallCancer.OnePolyp{x1, x2});
            R_Alive.OnePolyp{x1}              = cat(2, R_Alive.OnePolyp{x1}, Alive.OnePolyp{x1, x2});
            R_IntervallCancer.TwoPolyp{x1}    = cat(2, R_IntervallCancer.TwoPolyp{x1}, IntervallCancer.TwoPolyp{x1, x2});
            R_Alive.TwoPolyp{x1}              = cat(2, R_Alive.TwoPolyp{x1}, Alive.TwoPolyp{x1, x2});
            R_IntervallCancer.OneAdvPolyp{x1} = cat(2, R_IntervallCancer.OneAdvPolyp{x1}, IntervallCancer.OneAdvPolyp{x1, x2});
            R_Alive.OneAdvPolyp{x1}           = cat(2, R_Alive.OneAdvPolyp{x1}, Alive.OneAdvPolyp{x1, x2});
        end
    end
end

% we calculate mean and standard deviation
Readouts = {'All', 'NoPolyp', 'OnePolyp', 'TwoPolyp', 'OneAdvPolyp'};

for x1=15:length(Headlines) % we calculate an average
    for x2=1:length(Readouts)
        Mean_IntervallCancer.(Readouts{x2}){x1} = mean(R_IntervallCancer.(Readouts{x2}){x1});
        Mean_Alive.(Readouts{x2}){x1}           = mean(R_Alive.(Readouts{x2}){x1});
        
        Std_IntervallCancer.(Readouts{x2}){x1}  = std(R_IntervallCancer.(Readouts{x2}){x1});
        Std_Alive.(Readouts{x2}){x1}            = std(R_Alive.(Readouts{x2}){x1});
    end
end

% we write the Excel files
for x1=1:length(Readouts)
    xlswrite (ExcelFileName, Headlines, [Readouts{x1} '_ca_mean'], 'B1')
    xlswrite (ExcelFileName, Headlines, [Readouts{x1} '_ca_std'], 'B1')
    xlswrite (ExcelFileName, Headlines, [Readouts{x1} '_alive_mean'], 'B1')
    
    xlswrite (ExcelFileName, SummaryLegend, [Readouts{x1} '_ca_mean'], 'A2')
    xlswrite (ExcelFileName, SummaryLegend, [Readouts{x1} '_ca_std'], 'A2')
    xlswrite (ExcelFileName, SummaryLegend, [Readouts{x1} '_alive_mean'], 'A2')
    for f=15:length(Headlines)
        xlswrite (ExcelFileName, Mean_IntervallCancer.(Readouts{x1}){f}, [Readouts{x1} '_ca_mean'], HeadLineColumns{f})
        xlswrite (ExcelFileName, Std_IntervallCancer.(Readouts{x1}){f},  [Readouts{x1} '_ca_std'], HeadLineColumns{f})
        xlswrite (ExcelFileName, Mean_Alive.(Readouts{x1}){f},           [Readouts{x1} '_alive_mean'], HeadLineColumns{f})
    end
end

