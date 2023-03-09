function[elevation, reg, ice, soot, X, Y] = read_surface_data(output_dir_path)
%
% This function receives the path to the output directory of REGOLIT
% and returns the surface elevation, regolith fraction, ice fraction
% and "soot" fraction as 3-D NxNxM arrays, where N is the grid size 
% and M is the number of time steps.
%
% Written by Lior Rubanenko
% Technion Israel Institute of Science 2023
%

if nargin == 0
    output_dir_path = '/Users/liorr/Work/Projects/surface_evolution_model/REGOLIT/output/';
end

fidx = fopen([output_dir_path, 'x.out'], 'rb');
fidy = fopen([output_dir_path, 'y.out'], 'rb');

x = fread(fidx, [1, Inf], 'double');
y = fread(fidy, [1, Inf], 'double');
sz = [length(x), length(y)];

[X,Y] = meshgrid(x,y);

elevation_file_list = dir(fullfile(output_dir_path, 'elevation*.out'));
regolith_file_list = dir(fullfile(output_dir_path, 'regolith*.out'));
ice_file_list = dir(fullfile(output_dir_path, 'ice*.out'));
soot_file_list = dir(fullfile(output_dir_path, 'soot*.out'));

for ii = 1:length(elevation_file_list)
    elev_file_name = elevation_file_list(ii).name;
    reg_file_name = regolith_file_list(ii).name;
    ice_file_name = ice_file_list(ii).name;
    soot_file_name = soot_file_list(ii).name;
    
    elev_fid = fopen(fullfile(output_dir_path, elev_file_name), 'rb');
    reg_fid = fopen(fullfile(output_dir_path, reg_file_name), 'rb');
    ice_fid = fopen(fullfile(output_dir_path, ice_file_name), 'rb');
    soot_fid = fopen(fullfile(output_dir_path, soot_file_name), 'rb');
    
    elevation(:,:,ii) = fread(elev_fid, sz, 'double');
    reg(:,:,ii) = fread(reg_fid, sz, 'double');
    ice(:,:,ii) = fread(ice_fid, sz, 'double');
    soot(:,:,ii) = fread(soot_fid, sz, 'double');
    
    fclose(elev_fid);
    fclose(reg_fid);
    fclose(ice_fid);
    fclose(soot_fid);
end

fclose(fidx);
fclose(fidy);

return 