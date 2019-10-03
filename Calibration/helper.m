
N = length(initialCalibration.NewPolypRate);
sigmoidal = @(a0,a1,a2,x)(a0./(1+exp(-(a1.*x - a2))));
gaussian  = @(b0,b1,b2,x)(b0*exp(-(b1.*x - b2).^2));


%% dealing with new polyp rate
F = @(x)(sigmoidal(x(1),x(2),x(3),1:N) - initialCalibration.NewPolypRate);
xF = lsqnonlin(F,[0.03,0.1,1],[0 0 0],[]);


figure(1)
clf
hold on
plot(initialCalibration.NewPolypRate,'b')
plot(sigmoidal(xF(1),xF(2),xF(3),1:N),'r--')
hold off
xF

%% dealing with earlu progression rate
N = length(initialCalibration.EarlyProgressionRate);
F = @(x)(gaussian(x(1),x(2),x(3),1:N) - initialCalibration.EarlyProgressionRate);
xF = lsqnonlin(F,[1,0.1,0.6],[0 0 0],[]);

figure(2)
clf
hold on
plot(initialCalibration.EarlyProgressionRate,'b')
plot(gaussian(xF(1),xF(2),xF(3),1:N),'r--')
hold off
xF

%% dealing with advanced progression rate
N = length(initialCalibration.AdvancedProgressionRate);
F = @(x)(gaussian(x(1),x(2),x(3),1:N) - initialCalibration.AdvancedProgressionRate);
xF = lsqnonlin(F,[0.02,0.1,2.5],[0 0 0],[]);

figure(3)
clf
hold on
plot(initialCalibration.AdvancedProgressionRate,'b')
plot(gaussian(xF(1),xF(2),xF(3),1:N),'r--')
hold off
xF

%% checking individual risk settings
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


figure(1)
clf
hold on
plot(initialCalibration.IndividualRisk,'r')
plot(tmp,'b--')
hold off