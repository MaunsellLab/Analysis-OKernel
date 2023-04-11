
% JJC 20230411
% Use this to make patches for figures.
% calls gaborFig function

%   Parameters (Default):
%       theta (2*pi*rand): orientation of the gabor patch
%       lambda (20): spatial wavelength
%       Sigma (10): standard deviation of gaussian window
%       width (256): width of generated image
%       height (256): height of generated image
%       px (rand): horizontal center of gabor patch, relative to the image width (must be between 0 and 1)
%       py (rand): vertical center of gabor patch, relative to the image height (must be between 0 and 1)
%

sig = 10;
columns = 200;
rows = 200;

% Make Gabor Patch
[x, y, ~, F_odd] = gaborFig('theta', 0, 'lambda', 20,...
    'Sigma', sig, 'width', columns, 'height', rows, 'px', 0.5, 'py', 0.5);
pcolor(x,y,F_odd); axis image;
shading('interp'); colormap gray;

% Make Luminance Patch
[x, y, F_even, ~] = gaborFig('theta', 0, 'lambda', 1,...
    'Sigma', sig, 'width', columns, 'height', rows, 'px', 0.5, 'py', 0.5);
pcolor(x,y,-F_even); axis image;
shading('interp'); 
colormap gray;


