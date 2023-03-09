function[] = make_animated_gif(X, Y, data_array, file_path, cmap, coloraxis)
%
% This function receives a REGOLIT data arrays created using the read
% surface / read integrated surface functions, creates an animated gif
% from the data array and saves it in the given file_path
%
% Written by Lior Rubanenko, Technion Israel Institute of Technology
%
%

if nargin < 4
    error('Not enough input arguments');

elseif nargin == 4
    cmap = 'gray';
    coloraxis = [-1, 1];
    
elseif nargin == 5
    coloraxis = [-1, 1];
end

sz = size(data_array);

for ii = 1:sz(3)
    pcolor(X, Y, data_array(:,:,ii));
    colormap(cmap);
    axis equal
    axis tight
    shading flat
    caxis(coloraxis);
    
    frame = getframe(gca);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    
    if ii == 1
        imwrite(imind,cm,file_path,'gif', 'Loopcount',inf);
    else
        imwrite(imind,cm,file_path,'gif','WriteMode','append');
    end
    
end

end