function[params] = parse_config(config_filepath)
if nargin == 0
    config_filepath = '../config/config.cfg';
end

% Open the file and read the contents into a cell array
fid = fopen(config_filepath);
contents = textscan(fid, '%s', 'Delimiter', '\n');
fclose(fid);
contents = contents{1};

% Parse the file contents and extract the parameters
params = struct();
for i = 1:length(contents)
    line = strtrim(contents{i});
    
    % Skip comment lines
    if startsWith(line, '//')
        continue;
    end
    
    % Split the line into a variable name and a value
    [name, value] = strtok(line);
    
    if length(name) == 0
        continue;
    end
    
    value = str2double(strtrim(value));
    
    % Add the variable to the params struct
    params.(name) = value;
end
end