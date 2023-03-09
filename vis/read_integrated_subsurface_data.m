function[reg, ice, soot, X, Y] = read_integrated_subsurface_data(output_dir_path)
%
% This function receives the path to the output directory of REGOLIT
% and returns the integrated subsurface regolith fraction, ice fraction
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

regolith_file_list = dir(fullfile(output_dir_path, 'depthRegolith*.out'));
ice_file_list = dir(fullfile(output_dir_path, 'depthIce*.out'));
soot_file_list = dir(fullfile(output_dir_path, 'depthSoot*.out'));

for ii = 1:length(regolith_file_list)
    reg_file_name = regolith_file_list(ii).name;
    ice_file_name = ice_file_list(ii).name;
    soot_file_name = soot_file_list(ii).name;
    
    reg_fid = fopen(fullfile(output_dir_path, reg_file_name), 'rb');
    ice_fid = fopen(fullfile(output_dir_path, ice_file_name), 'rb');
    soot_fid = fopen(fullfile(output_dir_path, soot_file_name), 'rb');
    
    reg(:,:,ii) = fread(reg_fid, sz, 'double');
    ice(:,:,ii) = fread(ice_fid, sz, 'double');
    soot(:,:,ii) = fread(soot_fid, sz, 'double');
    
    fclose(reg_fid);
    fclose(ice_fid);
    fclose(soot_fid);
end

fclose(fidx);
fclose(fidy);

return 