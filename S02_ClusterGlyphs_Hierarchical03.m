

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
% vFilenames = { 'Rosetta Stone_All_Test.tif' };
vFilenames = { 'Rosetta Stone_All.tif' };
% vFilenames = { 'Rosetta Stone_All.tif' 'Roman_Font_Test.tif' };


mGlyphPolygonsAll = [];

for k = 1:length(vFilenames)

	strFilename = vFilenames{k}
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
% 	figure(1);
% 	clf;
% 	imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]));
% 	SaveNiceFigure;
% 	
% 	figure(2);
% 	clf;
% 	imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]));
% 	hold on;
% 	
% 	figure(3);
% 	clf;
% 	imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]));
% 	hold on;

	iPts = 100;

	mGlyphPolygons = zeros(iPts*2,nGlyphs);
	for i = 1:length(vGlyphs)
		vGlyphPolygon = vGlyphs{i};
		%plot(iWidth-vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 2);

		% Because the image was flipped, the polygon x values will be wrong
		vGlyphPolygon(:,1) = iWidth-vGlyphPolygon(:,1);
		
% 		figure(2); % Plot before interpolation
% 		plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), '.-w', 'LineWidth', 1, 'MarkerSize', 15);
		
		
		% Interpolate the polygon so that all glyphs have the same # of points
		vGlyphPolygon = interppoly(vGlyphPolygon, iPts);


% 		figure(3); % Plot after interpolation
		%Use this plot function when the image is transposed
% 		plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), '.-w', 'LineWidth', 1, 'MarkerSize', 15);
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
	
% 	figure(2);
% 	SaveNiceFigure;
% 	
% 	figure(3);
% 	SaveNiceFigure;
	
	mGlyphPolygonsAll = [ mGlyphPolygonsAll ; mGlyphPolygons ];
end

mGlyphPolygons = mGlyphPolygonsAll;

	
%% Get Principal Components and plot a subset of axes
[mTPC,mGPPC,~] = getPC(mGlyphPolygons);
%mTPC = eye(size(mTPC,1));
mGPPC = mGlyphPolygons * mTPC;


%% Use some subset of PCs to create a clustering tree

vDistMetrics = {
	'euclidean'					% 1
	'squaredeuclidean'			% 2
	'seuclidean'				% 3
	'cityblock'					% 4
	'minkowski'					% 5
	'chebychev'					% 6
	'mahalanobis'				% 7
	'cosine'					% 8
	'correlation'				% 9
	'spearman'					% 10
	'hamming'					% 11
	'jaccard'					% 12
};

j = [ 6 ];
nClusters = 500

X = mGPPC(:,1:20);

Y = pdist(X,vDistMetrics{j});
Z = linkage(Y);

T = cluster(Z,'maxclust',nClusters);


%%% Plot clusters
close all

% Sort clusters according to number of members and display the biggest
% first
nT = histc(T,1:nClusters);
[nT,iSortedClusters] = sort(nT,'descend');

plot(nT)
iStart = 1; %find(nT < 20, 1, 'first');
iEnd = find(nT < 4, 1, 'first');
for iKeep = iStart:iEnd
	i = iSortedClusters(iKeep);
	n = sum(T==i)
	
	figure(iKeep+100);
	
	vCluster = find(T==i);
	nI = length(vCluster);
	vCluster = vCluster(randperm(nI, min(nI,15)));
	
	ShowGlyph(mGlyphPolygons, vCluster);
% 		title(i);
	axis equal
	axis off
	set(gcf,'color','k');
	getframe;

	SaveNiceFigure;
	
	%[ i n ]
end

sprintf('%1.4f - %d: %s', skewness(T), j, vDistMetrics{j} )
%pause

%% Show first PCs in clusters
figure(2);
clf;
hold on;

vColors = iris( nClusters );
iStart = 2;%find(nT < 100, 1, 'first')
iEnd = nClusters%find(nT < 5, 1, 'first')
for i = iSortedClusters(iStart:iEnd)'
	n = sum(T==i)
	if( n > 1 )
		plot3d( mGPPC(find(T==i),1), mGPPC(find(T==i),2), mGPPC(find(T==i),3),...
			'.', 'Color', vColors(i,:), 'MarkerSize', 40 );
	end
	%[ i n ]
end

%%




