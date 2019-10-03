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

function [data, BM, err] = Evaluation_for_Jan_Sept_2019_lean(data, Variables, benchmarks, toPlot)

% comment about TIME (year)
% in all scripts PRECEDING this time was 1-100 (Lebensjahre)
% in EVALUATION we transform this to age (0-99)

% for PBC we save these data directly to the results variable
% Results = struct;
% Results.DiagnosedCancer = data.DiagnosedCancer;
% Results.DeathYear       = data.DeathYear;

%y = data.y;

% TumorRecordOut.Stage = zeros(10000, yearsToSimulate);
% TumorRecordOut.Location = zeros(10000, yearsToSimulate);
% TumorRecordOut.Sojourn = zeros(10000, yearsToSimulate);
% TumorRecordOut.DwellTime = zeros(10000, yearsToSimulate);
% TumorRecordOut.Gender = zeros(10000, yearsToSimulate);
% TumorRecordOut.Detection = zeros(10000, yearsToSimulate);
% TumorRecordOut.PatientNumber = zeros(10000, yearsToSimulate);
% TumorRecordOut.Time = zeros(10000, yearsToSimulate); % FIX BM 27.10.2018 
% 
% for i = 1:yearsToSimulate
%     idxTmp = (TumorRecord.Year == i);
%     if any(idxTmp)
%         Nrec = sum(idxTmp);
%         TumorRecordOut.Stage(1:Nrec,i) = TumorRecord.Stage(idxTmp);
%         TumorRecordOut.Location(1:Nrec,i) = TumorRecord.Location(idxTmp);
%         TumorRecordOut.Sojourn(1:Nrec,i) = double(TumorRecord.Sojourn(idxTmp))/double(numPeriods);
%         TumorRecordOut.DwellTime(1:Nrec,i) = double(TumorRecord.DwellTime(idxTmp))/double(numPeriods);
%         TumorRecordOut.Gender(1:Nrec,i) = TumorRecord.Gender(idxTmp);
%         TumorRecordOut.Detection(1:Nrec,i) = TumorRecord.Detection(idxTmp);
%         TumorRecordOut.PatientNumber(1:Nrec,i) = TumorRecord.SubjectID(idxTmp);
%         TumorRecordOut.Time(1:Nrec,i) = TumorRecord.Time(idxTmp); % FIX BM 27.10.2018 
%     end
% end


% flds = fields(TumorRecordOut);
% for i = 1:length(flds)
%     TumorRecordOut.(flds{i}) = TumorRecordOut.(flds{i})';
% end


%Variables.Benchmarks.Cancer.SymptomaticStageDistribution  = [15 35.6 27.9 21.5];
% now I just averaged the benchmarks of the 1988-2000 period
%Variables.Benchmarks.Cancer.SymptomaticStageDistribution  = %[18.92 27.67 29.89 23.52];

%Variables.Benchmarks.Cancer.ScreeningStageDistribution    = %[39.5 34.7 17.3 8.5];

%Variables.Benchmarks.Cancer.LocationRectumMale   = %[41.2     34.1      28.6     23.8];
%Variables.Benchmarks.Cancer.LocationRectumFemale = %[37.2     28.3      23.0     19.0];
%Variables.Benchmarks.Cancer.LocationRectumYear   = %{[51 55], [61 65], [71 75], [81 85]};  % year adapted

%Variables.Benchmarks.Cancer.Fastcancer           = %[0.005 0.05 0.08 0.25 3 20];

%SEER 1988 - 2000
%Variables.Benchmarks.Cancer.Ov_y_mort  =  %[1.5  5.5  12   17     22   27   32  37   42  47   52   57   62    67   72    77     82     87];
%Variables.Benchmarks.Cancer.Ov_mort     = %[0   0    0    0.1    0.2  0.4  1   2    4   8.2  15.8 28.6 47.2  71.2 102.6 140.4  193.5  285.1];
%Variables.Benchmarks.Cancer.Male_mort   = %[0   0    0    0.1    0.2  0.5  1.1 2.1  4.3 9.2  18.2 34   58.1  89.2 128.8 176.2  241.5  342.2];
%Variables.Benchmarks.Cancer.Female_mort = %[0   0    0    0.1    0.2  0.4  0.9 1.8  3.7 7.3  13.5 23.6 37.5  56.2 82.6  116.5  166.8  262.7];



% we calculate the number of patients with polpys 1-4 and express 
% them as percentage of survivors
% for f=1:y
%     NumPolyps(f)    = sum(data.MaxPolyps(f, data.YearIncluded(f, :)==1)>0); 
%     NumPolyps_2(f)  = sum(data.MaxPolyps(f, data.YearIncluded(f, :)==1)>1);
%     NumPolyps_3(f)  = sum(data.MaxPolyps(f, data.YearIncluded(f, :)==1)>2);
%     NumPolyps_4(f)  = sum(data.MaxPolyps(f, data.YearIncluded(f, :)==1)>3);
%     NumPolyps_5(f)  = sum(data.MaxPolyps(f, data.YearIncluded(f, :)==1)>4);
%     NumPolyps_6(f)  = sum(data.MaxPolyps(f, data.YearIncluded(f, :)==1)>5); 
% end
% for f=1:y
%     FracPolyps(f) = NumPolyps(f)/sum(data.YearIncluded(f, :))*100; 
%     FracPolyps_2(f) = NumPolyps_2(f)/sum(data.YearIncluded(f, :))*100;
%     FracPolyps_3(f) = NumPolyps_3(f)/sum(data.YearIncluded(f, :))*100;
%     FracPolyps_4(f) = NumPolyps_4(f)/sum(data.YearIncluded(f, :))*100;
%     FracPolyps_5(f) = NumPolyps_5(f)/sum(data.YearIncluded(f, :))*100;
%     FracPolyps_6(f) = NumPolyps_6(f)/sum(data.YearIncluded(f, :))*100;
% end

% the fraction of surviving patients with advanced polyps
Variables.Benchmarks.EarlyPolyp.Ov_y = benchmarks.EarlyPolyps.Ov_y';
Variables.Benchmarks.EarlyPolyp.Ov_perc = benchmarks.EarlyPolyps.Ov_perc';
FracPolyps = data.PolypsSumm.FracPolyps(1,:);
OutputValues = CalculateAgreement(FracPolyps, Variables.Benchmarks, 'EarlyPolyp', 'Ov_y', 'Ov_perc','Polyp'); 
BM.OutputValues.EarlyAdenoma_Ov = OutputValues;
err = (BM.OutputValues.EarlyAdenoma_Ov - Variables.Benchmarks.EarlyPolyp.Ov_perc)/max(Variables.Benchmarks.EarlyPolyp.Ov_perc);
if toPlot
   figure(1)
   clf
   subplot(4,3,1)
   hold on
    N = length(Variables.Benchmarks.EarlyPolyp.Ov_perc);
    plot(1:N,Variables.Benchmarks.EarlyPolyp.Ov_perc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.EarlyAdenoma_Ov);
   hold off
   title('early polyps overall')
end

% the fraction of surviving patients with advanced polyps
Variables.Benchmarks.AdvPolyp.Ov_y = benchmarks.AdvPolyps.Ov_y';
Variables.Benchmarks.AdvPolyp.Ov_perc = benchmarks.AdvPolyps.Ov_perc';
FracPolyps_5 = data.PolypsSumm.FracPolyps(5,:);
OutputValues = CalculateAgreement(FracPolyps_5, Variables.Benchmarks, 'AdvPolyp', 'Ov_y', 'Ov_perc', 'Polyp');  
BM.OutputValues.AdvAdenoma_Ov = OutputValues;

err = [err (BM.OutputValues.AdvAdenoma_Ov - Variables.Benchmarks.AdvPolyp.Ov_perc)/max(Variables.Benchmarks.AdvPolyp.Ov_perc)];
if toPlot
   subplot(4,3,2)
   hold on
    N = length(Variables.Benchmarks.AdvPolyp.Ov_perc);
    plot(1:N,Variables.Benchmarks.AdvPolyp.Ov_perc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.AdvAdenoma_Ov);
   hold off
   title('advanced polyps overall')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Early/ advanced polyps Male/ Female    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear EarlyPolyps AdvPolyps Included
% % we calculate the presence of polyps (all polyps or Advanced polyps and
% % express as percent of survivors
% for f1=1:2 % 1=male, 2=female
%     Gender = data.Gender == f1;
%     for f=1:y
%         EarlyPolyps{f1}(f) = sum(data.MaxPolyps(f, and(Gender, data.YearIncluded(f, :)))>0); 
%         AdvPolyps{f1}(f)   = sum(data.MaxPolyps(f, and(Gender, data.YearIncluded(f, :)))>4);
%         Included           = sum(data.YearIncluded(f, Gender));
%         EarlyPolyps{f1}(f) = EarlyPolyps{f1}(f)/Included*100; 
%         AdvPolyps{f1}(f)   = AdvPolyps{f1}(f)/Included*100; 
%     end
% end

% Early polyps male
Variables.Benchmarks.EarlyPolyp.Male_y = benchmarks.EarlyPolyps.Male_y';
Variables.Benchmarks.EarlyPolyp.Male_perc = benchmarks.EarlyPolyps.Male_perc';
OutputValues = CalculateAgreement(data.PolypsSumm.FracPolypsMale(1,:), Variables.Benchmarks, 'EarlyPolyp', 'Male_y', 'Male_perc','Polyp'); 
BM.OutputValues.EarlyAdenoma_Male = OutputValues;
err = [err (BM.OutputValues.EarlyAdenoma_Male - Variables.Benchmarks.EarlyPolyp.Male_perc)/max(Variables.Benchmarks.EarlyPolyp.Male_perc)];

if toPlot
   subplot(4,3,3)
   hold on
    N = length(Variables.Benchmarks.EarlyPolyp.Male_perc);
    plot(1:N,Variables.Benchmarks.EarlyPolyp.Male_perc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.EarlyAdenoma_Male);
   hold off
   title('early polyps male')
end

% Early polyps female
Variables.Benchmarks.EarlyPolyp.Female_y = benchmarks.EarlyPolyps.Female_y';
Variables.Benchmarks.EarlyPolyp.Female_perc = benchmarks.EarlyPolyps.Female_perc';
OutputValues = CalculateAgreement(data.PolypsSumm.FracPolypsFemale(1,:), Variables.Benchmarks, 'EarlyPolyp', 'Female_y', 'Female_perc', 'Polyp'); 
BM.OutputValues.EarlyAdenoma_Female = OutputValues;
err = [err (BM.OutputValues.EarlyAdenoma_Female - Variables.Benchmarks.EarlyPolyp.Female_perc)/Variables.Benchmarks.EarlyPolyp.Female_perc];
if toPlot
   subplot(4,3,4)
   hold on
    N = length(Variables.Benchmarks.EarlyPolyp.Female_perc);
    plot(1:N,Variables.Benchmarks.EarlyPolyp.Female_perc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.EarlyAdenoma_Female);
   hold off
   title('early polyps female')
end

% advanced polyps male
Variables.Benchmarks.AdvPolyp.Male_y = benchmarks.AdvPolyps.Male_y';
Variables.Benchmarks.AdvPolyp.Male_perc = benchmarks.AdvPolyps.Male_perc';
OutputValues = CalculateAgreement(data.PolypsSumm.FracPolypsMale(5,:), Variables.Benchmarks, 'AdvPolyp', 'Male_y', 'Male_perc', 'Polyp'); 
BM.OutputValues.AdvAdenoma_Male = OutputValues;
err = [err (BM.OutputValues.AdvAdenoma_Male - Variables.Benchmarks.AdvPolyp.Male_perc)/max(Variables.Benchmarks.AdvPolyp.Male_perc)];
if toPlot
   subplot(4,3,5)
   hold on
    N = length(Variables.Benchmarks.AdvPolyp.Male_perc);
    plot(1:N,Variables.Benchmarks.AdvPolyp.Male_perc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.AdvAdenoma_Male);
   hold off
   title('advanced polyps male')
end

% advanced polyps female
Variables.Benchmarks.AdvPolyp.Female_y = benchmarks.AdvPolyps.Female_y';
Variables.Benchmarks.AdvPolyp.Female_perc = benchmarks.AdvPolyps.Female_perc';
OutputValues = CalculateAgreement(data.PolypsSumm.FracPolypsFemale(5,:), Variables.Benchmarks, 'AdvPolyp', 'Female_y', 'Female_perc', 'Polyp'); 
BM.OutputValues.AdvAdenoma_Female = OutputValues;
err = [err (BM.OutputValues.AdvAdenoma_Female - Variables.Benchmarks.AdvPolyp.Female_perc)/max(Variables.Benchmarks.AdvPolyp.Female_perc)];
if toPlot
   subplot(4,3,6)
   hold on
    N = length(Variables.Benchmarks.AdvPolyp.Female_perc);
    plot(1:N,Variables.Benchmarks.AdvPolyp.Female_perc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.AdvAdenoma_Female);
   hold off
   title('advanced polyps female')
end
clear AdvPolyps EarlyPolyps Gender Included


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Cancer Incidence All  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SumIncidenceCounter = [sum(data.IncidenceCounter(:,1:4),2) sum(data.IncidenceCounter(:,5:8),2)   sum(data.IncidenceCounter(:,11:15),2) ...
    sum(data.IncidenceCounter(:,16:20),2) sum(data.IncidenceCounter(:,21:25),2) sum(data.IncidenceCounter(:,26:30),2)... % year adapted
       sum(data.IncidenceCounter(:,31:35),2)   sum(data.IncidenceCounter(:,36:40),2) sum(data.IncidenceCounter(:,41:45),2) ...
       sum(data.IncidenceCounter(:,46:50),2) sum(data.IncidenceCounter(:,51:55),2) sum(data.IncidenceCounter(:,56:60),2)...
       sum(data.IncidenceCounter(:,61:65),2)   sum(data.IncidenceCounter(:,66:70),2) sum(data.IncidenceCounter(:,71:75),2) ...
       sum(data.IncidenceCounter(:,76:80),2) sum(data.IncidenceCounter(:,81:85),2) sum(data.IncidenceCounter(:,86:90),2)];

%for f=1:y
%    i(f) = length(find(data.TumorRecord.Stage(f, :)));
%    j(f) = sum(data.YearIncluded(f, :));
%end  
% we summarize in 5 year intervals
% SumCa =  [sum(i(1:4)) sum(i(5:8))   sum(i(11:15))  sum(i(16:20)) sum(i(21:25)) sum(i(26:30))... % year adapted
%       sum(i(31:35))   sum(i(36:40)) sum(i(41:45)) sum(i(46:50)) sum(i(51:55)) sum(i(56:60))...
%       sum(i(61:65))   sum(i(66:70)) sum(i(71:75)) sum(i(76:80)) sum(i(81:85)) sum(i(86:90))];
% SumPat = [sum(j(1:4)) sum(j(5:8))   sum(j(11:15))  sum(j(16:20)) sum(j(21:25)) sum(j(26:30))...
%       sum(j(31:35))   sum(j(36:40)) sum(j(41:45)) sum(j(46:50)) sum(j(51:55)) sum(j(56:60))...
%       sum(j(61:65))   sum(j(66:70)) sum(j(71:75)) sum(j(76:80)) sum(j(81:85)) sum(j(86:90))];
%   
% % and express is as new cancer cases per 100'000 patients
% for f=1:length(SumCa)
%     Incidence(f) = SumCa(f)/SumPat(f);
% end
% Incidence = Incidence * 100000;


% Overall cancer incidence
Variables.Benchmarks.Cancer.Ov_y = benchmarks.Cancer.Ov_y';
Variables.Benchmarks.Cancer.Ov_inc = benchmarks.Cancer.Ov_inc';
OutputValues = CalculateAgreement(SumIncidenceCounter(1,:)./SumIncidenceCounter(4,:)*100000, Variables.Benchmarks, 'Cancer', 'Ov_y', 'Ov_inc', 'Cancer');  
BM.OutputValues.Cancer_Ov = OutputValues;
err = [err (BM.OutputValues.Cancer_Ov  - Variables.Benchmarks.Cancer.Ov_inc)/max(Variables.Benchmarks.Cancer.Ov_inc)];
if toPlot
   subplot(4,3,7)
   hold on
    N = length(Variables.Benchmarks.Cancer.Ov_inc);
    plot(1:N,Variables.Benchmarks.Cancer.Ov_inc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.Cancer_Ov);
   hold off
   title('overall cancer incidence')
end

% BM.Incidence = Incidence; 
clear Incidence SumCa SumPat i j 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   Cancer Incidence Male/ Female   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear Counter Incidence1 i j tmp3 tmp4 tmp5
 
% for f1=1:2
%     for f2=1:y
%         i(f2) = length(find(data.TumorRecord.Stage(f2, data.TumorRecord.Gender(f2, :)==f1)));
%         j(f2) = sum(data.YearIncluded(f2, data.Gender ==f1));
%     end    
%     tmp3 =  [sum(i(1:4)) sum(i(5:8))   sum(i(11:15))  sum(i(16:20)) sum(i(21:25)) sum(i(26:30))... % year adapted
%       sum(i(31:35))   sum(i(36:40)) sum(i(41:45)) sum(i(46:50)) sum(i(51:55)) sum(i(56:60))...
%       sum(i(61:65))   sum(i(66:70)) sum(i(71:75)) sum(i(76:80)) sum(i(81:85)) sum(i(86:90))];
%     tmp4 = [sum(j(1:4)) sum(j(5:8))   sum(j(11:15))  sum(j(16:20)) sum(j(21:25)) sum(j(26:30))...
%       sum(j(31:35))   sum(j(36:40)) sum(j(41:45)) sum(j(46:50)) sum(j(51:55)) sum(j(56:60))...
%       sum(j(61:65))   sum(j(66:70)) sum(j(71:75)) sum(j(76:80)) sum(j(81:85)) sum(j(86:90))];
%   
%     for f=1:length(tmp3)
%         tmp5(f) = tmp3(f)/tmp4(f);
%     end
%     Incidence{f1} = tmp5 *100000;
% end

% male cancer incidence

Variables.Benchmarks.Cancer.Male_y = benchmarks.Cancer.Male_y';
Variables.Benchmarks.Cancer.Male_inc = benchmarks.Cancer.Male_inc';
OutputValues = CalculateAgreement(SumIncidenceCounter(2,:)./SumIncidenceCounter(5,:)*100000, Variables.Benchmarks, 'Cancer', 'Male_y', 'Male_inc', 'Cancer');  
BM.OutputValues.Cancer_Male = OutputValues;
err = [err (BM.OutputValues.Cancer_Male  - Variables.Benchmarks.Cancer.Male_inc)/max(Variables.Benchmarks.Cancer.Male_inc)];
if toPlot
   subplot(4,3,8)
   hold on
    N = length(Variables.Benchmarks.Cancer.Male_inc);
    plot(1:N,Variables.Benchmarks.Cancer.Male_inc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.Cancer_Male);
   hold off
   title('male cancer incidence')
end

% female cancer incidence
Variables.Benchmarks.Cancer.Female_y = benchmarks.Cancer.Female_y';
Variables.Benchmarks.Cancer.Female_inc = benchmarks.Cancer.Female_inc';
OutputValues = CalculateAgreement(SumIncidenceCounter(3,:)./SumIncidenceCounter(6,:)*100000, Variables.Benchmarks, 'Cancer', 'Female_y', 'Female_inc', 'Cancer');  
BM.OutputValues.Cancer_Female = OutputValues;
err = [err (BM.OutputValues.Cancer_Female  - Variables.Benchmarks.Cancer.Female_inc)/max(Variables.Benchmarks.Cancer.Female_inc)];
if toPlot
   subplot(4,3,9)
   hold on
    N = length(Variables.Benchmarks.Cancer.Female_inc);
    plot(1:N,Variables.Benchmarks.Cancer.Female_inc,'Marker','o','LineStyle','none');
    plot(1:N, BM.OutputValues.Cancer_Female);
   hold off
   title('female cancer incidence')
end

clear Incidence i j tmp3 tmp4 tmp5 



% Polyp_early         = zeros(6,1);
% Polyp_adv           = zeros(6,1);
% BM_value_early      = zeros(6,1);
% BM_value_adv        = zeros(6,1);
% 
% BM_value            = Variables.Benchmarks.Polyp_Distr;
% 
% Summe_early         = sum(sum(data.AllPolyps(1:4, 51:76))); % year adapted
% Summe_adv           = sum(sum(data.AllPolyps(5:6, 51:76))); % year adapted
% BM_value_early(1:4) = BM_value(1:4)/sum(BM_value(1:4))*100;
% BM_value_adv(5:6)   = BM_value(5:6)/sum(BM_value(5:6))*100;

% for f=1:4
% %     BM.description{bmc} = ['% of all early polyps P ' num2str(f)];
%     Polyp_early(f) = sum(data.AllPolyps(f, 51:76))/Summe_early*100; % year adapted
% % %     BM.value{bmc} = Polyp_early(f); BM.benchmark{bmc} = BM_value_early(f);
% %     
% %     if isequal(f, 1), LinePos(f) = Polyp_early(f)/2;
% %     else LinePos(f) = sum(Polyp_early(1:f-1))+Polyp_early(f)/2;
% %     end
% %     if and(BM.value{bmc} > BM_value_early(f)*(1 - tolerance), BM.value{bmc} < (BM_value_early(f)*(1 + tolerance)))
% %         BM.flag{bmc} = 'green'; Color{f} = 'g';
% %     else
% %         BM.flag{bmc} = 'red';   Color{f} = 'r';
% %     end
% %     BM.Polyp_Distr(f) = BM.value{bmc};
% %     bmc = bmc +1;
% end
% for f=5:6
% %     BM.description{bmc} = ['% of all early polyps P ' num2str(f)];
%     Polyp_adv(f) = sum(data.AllPolyps(f, 51:76))/Summe_adv*100; % year adapted
% %     BM.value{bmc} = Polyp_adv(f); BM.benchmark{bmc} = BM_value_adv(f);
% %     
% %     if isequal(f, 1), LinePos(f) = Polyp_adv(f)/2;
% %     else LinePos(f) = sum(Polyp_adv(1:f-1))+Polyp_adv(f)/2;
% %     end
% %     if and(BM.value{bmc} > BM_value_adv(f)*(1 - tolerance), BM.value{bmc} < (BM_value_adv(f)*(1 + tolerance)))
% %         BM.flag{bmc} = 'green'; Color{f} = 'g';
% %     else
% %         BM.flag{bmc} = 'red';   Color{f} = 'r';
% %     end
% %     BM.Polyp_Distr(f) = BM.value{bmc};
% %     bmc = bmc +1;
% end
% % if DispFlag
% %     figure(h2), subplot(3,3,2)
% %     bar(cat(2, Polyp_early, zeros(6,1), BM_value_early, zeros(6,1), ...
% %         Polyp_adv, zeros(6,1), BM_value_adv)', 'stacked'), hold on
% %     for f=1:4, line([1.5 2.5], [LinePos(f) LinePos(f)], 'color', Color{f}), end
% %     for f=5:6, line([5.5 6.5], [LinePos(f) LinePos(f)], 'color', Color{f}), end
% %     l=legend('Adenoma 3mm', 'Adenoma 5mm', 'Adenoma 7mm', 'Adenoma 9mm', 'Adv Adenoma P5', 'Adv Adenoma P6');
% %     set(l, 'location', 'northoutside', 'fontsize', 6)
% %     ylabel('% of adenomas', 'fontsize', 6)
% %     set(gca, 'xticklabel', {'Ear.Ad.' '' 'BM' '' 'Adv.Ad.' '' 'BM'}, 'fontsize', 6, 'ylim', [0 100])
% % end
% BM.Polyp_early    = Polyp_early;
% % BM.BM_value_early = BM_value_early;
% BM.Polyp_adv      = Polyp_adv;
% % BM.BM_value_adv   = BM_value_adv;
% % BM.Pflag          = Color;
% clear LinePos Polyp Color Summe


% we summarize the number of polyps
% for f=1:y
%     FivePolyps(f) = sum(data.NumPolyps(f,:) > 4); 
%     FourPolyps(f) = sum(data.NumPolyps(f,:) > 3);
%     ThreePolyps(f)= sum(data.NumPolyps(f,:) > 2); 
%     TwoPolyps(f)  = sum(data.NumPolyps(f,:) > 1); 
%     OnePolyp(f)   = sum(data.NumPolyps(f,:) > 0); 
% end

% these data are for the next plot which uses uncorrected numbers (at least
% one polyp... we summarize the population of different ages
%NumYoung = 0; NumMid = 0; NumOld = 0; NumAllAges = 0;
%for f=41:55;  NumYoung = NumYoung+sum(data.YearIncluded(f, :)); end % year adapted
%for f=56:75;  NumMid   = NumMid+sum(data.YearIncluded(f, :)); end   % year adapted
%for f=76:91; NumOld   = NumOld+sum(data.YearIncluded(f, :)); end    % year adapted
%for f=50:100;  NumAllAges = NumAllAges+sum(data.YearIncluded(f, :)); end % year adapted

NumYoung = sum(data.IncidenceCounter(4,41:55));
NumMid = sum(data.IncidenceCounter(4,56:75));
NumOld = sum(data.IncidenceCounter(4,76:91));

YoungPop = transpose(sum(data.PolypsSumm.ActNumPolyps(:,41:55),2)./NumYoung*100);
MidPop = transpose(sum(data.PolypsSumm.ActNumPolyps(:,56:75),2)./NumMid*100);
OldPop = transpose(sum(data.PolypsSumm.ActNumPolyps(:,76:91),2)./NumOld*100);


BM.YoungPop=YoungPop;BM.MidPop=MidPop;BM.OldPop=OldPop;

% % we correct for multiple polyps
% AllPolyps   = OnePolyp(1:100)/100;
% OnePolyp    = OnePolyp - TwoPolyps;
% TwoPolyps   = TwoPolyps - ThreePolyps;
% ThreePolyps = ThreePolyps - FourPolyps;
% FourPolyps  = FourPolyps  - FivePolyps;

% if DispFlag
%     figure(h2), subplot(3,3,3)
%     plot(0:99, OnePolyp(1:100)./AllPolyps, 'color', 'r'), hold on
%     plot(0:99, TwoPolyps(1:100)./AllPolyps, 'color', 'k')
%     plot(0:99, ThreePolyps(1:100)./AllPolyps, 'color', 'b')
%     plot(0:99, FourPolyps(1:100)./AllPolyps, 'color', 'g')
%     plot(0:99, FivePolyps(1:100)./AllPolyps, 'color', 'm')
%     set(gca, 'xlim', [0 100],  'fontsize', FontSz)
%     xlabel('year'), ylabel('% of patients with adenomas', 'Fontsize', FontSz), title('number of adenomas', 'Fontsize', FontSz)
%     l=legend('1 adenoma', '2 adenomas', '3 adenomas', '4 adenomas', '>4 adenomas');
%     set(l, 'location', 'northoutside', 'fontsize', FontSz-1)
%     set(gca, 'xlim', [0 100],  'fontsize', FontSz)
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%    Number Polyps Frequency distribution    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Variables.Benchmarks.MultiplePolypsYoung = benchmarks.MultiplePolyps.MultiplePolypsYoung';
Variables.Benchmarks.MultiplePolypsOld   = benchmarks.MultiplePolyps.MultiplePolypsOld';
Variables.Benchmarks.MultiplePolyp = benchmarks.MultiplePolyps.MultiplePolyps';

YoungBenchmark = Variables.Benchmarks.MultiplePolypsYoung;
MidBenchmark   = Variables.Benchmarks.MultiplePolyp;
OldBenchmark   = Variables.Benchmarks.MultiplePolypsOld;


BM.OutputValues.YoungPop = YoungPop; BM.OutputValues.MidPop = MidPop; BM.OutputValues.OldPop = OldPop; 

err = [err (YoungBenchmark - YoungPop)/max(YoungBenchmark) ...
    (MidBenchmark - MidPop)/max(MidBenchmark) ....
    (OldBenchmark - OldPop)/max(OldBenchmark)];
if toPlot
   subplot(4,3,10)
   hold on
    N = length(YoungBenchmark);
    plot(1:N,YoungBenchmark,'Marker','o','LineStyle','none','Color','r');
        plot(1:N,MidBenchmark,'Marker','o','LineStyle','none','Color','b');
        plot(1:N,OldBenchmark,'Marker','o','LineStyle','none','Color','k');
    plot(1:N, YoungPop,'Color','r');
    plot(1:N, MidPop,'Color','b');
        plot(1:N, OldPop,'Color','k');
   hold off
   title('number polyps frequency distribution')
   legend({'young','mid','old'})
   
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%    Direct Cancer                           %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% for f=1:5
%     value(f) = sum(data.DirectCancer(f, 1:100))/sum(data.AllPolyps(f,1:100))*100;
% end
% value(6) = sum(data.ProgressedCancer(1:100))/sum(data.AllPolyps(6,1:100))*100;
% 
% % we correct to relative danger
% value = value./ sum(value)*100;
% BM.CancerOriginValue   = value;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%    Location           %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% tmp_all    = data.TumorRecord.Stage;
% tmp_male   = data.TumorRecord.Gender == 1;
% tmp_female = data.TumorRecord.Gender == 2;
% 
% tmp_Rectum = data.TumorRecord.Stage;
% tmp_Rectum(data.TumorRecord.Location <13) = 0;
% tmp_Right  = data.TumorRecord.Stage;
% tmp_Right(data.TumorRecord.Location >3) = 0;
% tmp_Rest   = data.TumorRecord.Stage;
% tmp_Rest(data.TumorRecord.Location ==13) = 0;
% tmp_Rest(data.TumorRecord.Location <4) = 0;
% 
% for f=1:4
%     Sum_Stage_all(f)    = sum(sum(tmp_all==f+6));
%     Sum_Stage_Rectum(f) = sum(sum(tmp_Rectum==f+6));
%     Sum_Stage_Right(f)  = sum(sum(tmp_Right==f+6));
%     Sum_Stage_Rest(f)   = sum(sum(tmp_Rest==f+6));
% end
% 
% tmp_Rectum_male  = tmp_Rectum>0;
% tmp_Rectum_male(tmp_female) = 0;
% tmp_Rest_male    = tmp_Rest>0;
% tmp_Rectum_male(tmp_female) = 0;
% tmp_all_male     = tmp_all >0;%m
% tmp_all_male(tmp_female) = 0; %m
% 
% tmp_Rectum_female  = tmp_Rectum>0;
% tmp_Rectum_female(tmp_male) = 0;
% tmp_Rest_female    = tmp_Rest>0;
% tmp_Rectum_female(tmp_male) = 0;
% tmp_all_female     = tmp_all >0;%m
% tmp_all_female(tmp_male) = 0; %m
% 
% clear value x 
% for f=1:100
%     LocationRectum{1}(f) = sum(tmp_Rectum_male(f,:));
%     LocationRest{1}(f)   = sum(tmp_Rest_male(f,:));
%     LocationRectum{2}(f) = sum(tmp_Rectum_female(f,:));
%     LocationRest{2}(f)   = sum(tmp_Rest_female(f,:));
%     LocationAll{1}(f)    = sum(tmp_all_male(f,:));
%     LocationAll{2}(f)    = sum(tmp_all_female(f,:));
% end
% 
% %%% benchmarks
% LocX               = Variables.Benchmarks.Cancer.LocationRectumYear;
%                     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% carcinoma rectum both genders                     %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % we average the male and female benchmarks
% % here we only collect the data for display during adjustment of adenomas
% LocationRectumAllGender = (LocationRectum{1}(1:100) + LocationRectum{2}(1:100))/2;
% LocationRest            = (LocationRest{1}(1:100) + LocationRest{2}(1:100))/2;
% 
% for f=1:length(LocX)
%     value(f) = sum(LocationRectumAllGender((LocX{f}(1)-2):(LocX{f}(2)+2)))/...
%         (sum(LocationRectumAllGender((LocX{f}(1)-2):(LocX{f}(2)+2))) + sum(LocationRest((LocX{f}(1)-2):(LocX{f}(2)+2))))*100;
%     BM.LocationRectum(f)     = value(f);
% end

%disp('ready')     

function OutputValues = CalculateAgreement(DataGraph, Benchmarks, Struct1, Struct2, Struct3, Flag) 
BM_year  = Benchmarks.(Struct1).(Struct2);
% BM_value = Benchmarks.(Struct1).(Struct3);

OutputValues = zeros(1, length(BM_year));

% add bench marks

for f=1:length(BM_year)
    if (BM_year(f) >5) && (BM_year(f) <95)
        if isequal(Flag, 'Polyp')
            value = mean(DataGraph(BM_year(f)-1 : BM_year(f)+3)); % year adapted
            OutputValues(f)= value; 
        elseif isequal(Flag, 'Cancer')
            if BM_year(f) >20  % we ignore benchmarks for age 1-20
                value = DataGraph(f); % year adapted
                OutputValues(f)= value; 
            end
        else
            error('wrong flag')
        end
    end
end