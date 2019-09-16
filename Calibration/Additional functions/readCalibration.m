function settings = readCalibration(file)
    if isempty(file)
       [file, path] = uigetfile({'*.mat'},...
                          'Select calibration file'); 
    end
    settings = load([path file]);
    settings = settings.Variables;
end

