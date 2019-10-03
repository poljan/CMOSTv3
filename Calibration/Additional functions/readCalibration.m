function [settings, path, file] = readCalibration(file)
    if isempty(file)
       [file, path] = uigetfile({'*.mat'},...
                          'Select calibration file'); 
    end
    settings = load([path file]);
    if isfield(settings,'Variables')
        settings = settings.Variables;
    else
       settings = settings.temp; 
    end
       
end

