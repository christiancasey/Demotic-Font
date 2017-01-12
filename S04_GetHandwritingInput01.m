

%% Script Initialization
switchToCD;
clc;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap jet;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure')

%% Initialize global variables and setup callbacks

figure;
plot([-1 101],[-1 101],'.');
axis([0 100 0 100]);
axis off;
P = get_pencil_curve;

img = frame2im( getframe(gcf) );
img = img(:,:,1);
img = (img > 128);
img = ~imdilate( ~img, strel( 'disk', 20 ) );
imshow(img);














