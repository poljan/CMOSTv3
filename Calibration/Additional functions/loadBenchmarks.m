function benchmarks = loadBenchmarks(folder)
    if isempty(folder) %prompt for folder with benchmarks
        folder = uigetdir([],'Select forlder with benchmarks');
    end
    files = dir([folder '/*.xlsx']);
    for i = 1:length(files)
        benchmarks.(files(i).name(1:end-5)) = readtable([folder '/' files(i).name]);
    end
end

