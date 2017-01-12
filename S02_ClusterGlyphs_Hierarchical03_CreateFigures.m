

%% Script Initialization
switchToCD;
clc;
clc;	clear;	close all;	
figure; figure(gcf);	whitebg('w');	colormap jet;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure')

%% Load image and begin processing

% vFilenames = { 'Roman_Font_Test_Small.tif' };
vFilenames = { 'Roman_Font_Test.tif' };
vFilenames = { 'Rosetta Stone_All_Test.tif' };
% vFilenames = { 'Rosetta Stone_All.tif' };
% vFilenames = { 'Rosetta Stone_All.tif' 'Roman_Font_Test.tif' };


mGlyphPolygonsAll = [];

% for k = 1:length(vFilenames)

	strFilename = vFilenames{1}
	img = imread(strFilename);
	[iHeight iWidth] = size(img);

	% Convert to logical and invert so that bwlabel works properly
	img = ~(img > 0);

	% Remove very small blobs, since these are probably just noise
	img = bwareaopen(img, 500);
	%imagescz(img);

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
	clf;
	imshow(label2rgb(imgGlyphs, @iris, [0 0 0]));
	set(gcf,'color','k');
	SaveNiceFigure;
	
	%%
	figure(2);
	clf;
	imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]*0));
	hold on;
	
	figure(3);
	clf;
	imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]*0));
	hold on;

	iPts = 11;

	mGlyphPolygons = zeros(iPts*2,nGlyphs);
	for i = 1:length(vGlyphs)
		vGlyphPolygon = vGlyphs{i};
		%plot(iWidth-vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 2);

		% Because the image was flipped, the polygon x values will be wrong
		vGlyphPolygon(:,1) = iWidth-vGlyphPolygon(:,1);
		
		figure(2); % Plot before interpolation
		plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), '.-w', 'LineWidth', 1, 'MarkerSize', 15);
		
		
		% Interpolate the polygon so that all glyphs have the same # of points
		vGlyphPolygon = interppoly(vGlyphPolygon, iPts);


		figure(3); % Plot after interpolation
		%Use this plot function when the image is transposed
		plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), '.-w', 'LineWidth', 1, 'MarkerSize', 15);
		% Use this plot to leave the original image unchanged
		%plot(vGlyphPolygon(:,2), vGlyphPolygon(:,1), 'w', 'LineWidth', 1);

		% Subtract the centroid to zero out polygons and allow analysis
		vGPCenter = mean(vGlyphPolygon);
		for j = 1:2
			vGlyphPolygon(:,j) = vGlyphPolygon(:,j)-vGPCenter(j);
		end

		%text(vGPCenter(1),vGPCenter(2),num2str(i),'Color','k','BackgroundColor','w');

		% Add the glyphs polygons to a big matrix
		mGlyphPolygons(:,i) = vGlyphPolygon(:);
	end

	% Transpose mGlyphPolygons so that they're oriented properly 
	% [Glyphs Points]
	mGlyphPolygons = mGlyphPolygons';
	
	figure(2);
	set(gcf,'color','k');
	SaveNiceFigure;
	
	figure(3);
	set(gcf,'color','k');
	SaveNiceFigure;
	
% 	mGlyphPolygonsAll = [ mGlyphPolygonsAll ; mGlyphPolygons ];
% end

% mGlyphPolygons = mGlyphPolygonsAll;


