

%% Script Initialization
switchToCD;
clc;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap jet;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure')

%% Load image and begin processing

strFilename = 'Roman_Font_Test_Small.tif';
strFilename = 'Roman_Font_Test.tif';
%strFilename = 'Rosetta Stone_All_Test.tif';
%strFilename = 'Rosetta Stone_All.tif';

img = imread(strFilename);
[iHeight iWidth] = size(img);

% Convert to logical and invert so that bwlabel works properly
img = ~(img > 0);

% Flip the image horizontally and then transpose so that labels are in
% approximately the right order. This is not a perfect solution.
img = fliplr(img)';

% Label blobs and get the number of distinct blobs for later looping
%[imgGlyphs nGlyphs] = bwlabel(img);

[vGlyphs,imgGlyphs] = bwboundaries(img,'noholes');
nGlyphs = max(imgGlyphs(:));

% Return both images to their orignal orientations.
img = fliplr(img');
imgGlyphs = fliplr(imgGlyphs');

%% Display the labeled image and polygons
figure(1);
imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]));
hold on;

iPts = 100;

mGlyphPolygons = zeros(iPts*2,nGlyphs);
for i = 1:length(vGlyphs)
	vGlyphPolygon = vGlyphs{i};
	%plot(iWidth-vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 2);
	
	% Interpolate the polygon so that all glyphs have the same # of points
	vGlyphPolygon = interppoly(vGlyphPolygon, iPts);
	
	% Because the image was flipped, the polygon x values will be wrong
	vGlyphPolygon(:,1) = iWidth-vGlyphPolygon(:,1);
	
	% Use this plot function when the image is transposed
	plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 1);
	% Use this plot to leave the original image unchanged
	%plot(vGlyphPolygon(:,2), vGlyphPolygon(:,1), 'w', 'LineWidth', 1);
	
	% Subtract the centroid to zero out polygons and allow analysis
	vGPCenter = mean(vGlyphPolygon);
	for j = 1:2
		vGlyphPolygon(:,j) = vGlyphPolygon(:,j)-vGPCenter(j);
	end
	
	text(vGPCenter(1),vGPCenter(2),num2str(i),'Color','k','BackgroundColor','w');
	
	% Add the glyphs polygons to a big matrix
	mGlyphPolygons(:,i) = vGlyphPolygon(:);
end

% Transpose mGlyphPolygons so that they're oriented properly
mGlyphPolygons = mGlyphPolygons';

title(iPts);
hold off;
axis normal;
axis equal;

	
%% Get Principal Components and plot a subset of axes
[~,mGPPC,~] = getPC(mGlyphPolygons);

%mGPPC = mGPPC';
iAxis1 = 1;
iAxis2 = 2;
iAxis3 = 3;

figure(21);
clf;
hold on;
plotxyz(mGPPC(:,[iAxis1 iAxis2 iAxis3]),'.w');
for i = 1:nGlyphs	
	text(mGPPC(i,iAxis1),mGPPC(i,iAxis2),mGPPC(i,iAxis3),num2str(i),'Color','k','BackgroundColor','w');
end
hold off;


%% Use some subset of PCs to create a clustering tree

X = mGPPC(:,1:10);
Y = pdist(X,'cityblock');
sqY = squareform(Y);



Z = linkage(Y);

figure(31);
dendrogram(Z,0);
zoom on;

%% Compare with the dendrogram created by the raw values
% 
% X = mGlyphPolygons;
% Y = pdist(X);
% %sqY = squareform(Y);
% 
% 
% 
% Z = linkage(Y);
% 
% figure(41);
% dendrogram(Z,0);
% zoom on;

%% Cluster

clc

% for c = 0.5:0.01:1.5

	T = cluster(Z,'maxclust',30)
% 	[ c max(T(:)) ]

% end

%% Plot clusters
close all

for i = 1:max(T)
	ShowGlyph(mGlyphPolygons, find(T==i));
	title(i);
	getframe;
end
