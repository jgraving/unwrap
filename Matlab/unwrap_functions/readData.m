function dt = readData(dataDir, file)
    % Read a dataset, with individual IDs as numbers
    dt = readtable(fullfile(dataDir, file));
    if ismember('ID', dt.Properties.VariableNames) && ~isnumeric(dt.ID)
        dt.ID = str2double(string(dt.ID));
    end
end
