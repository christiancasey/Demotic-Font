

%% Script Initialization
switchToCD;
clc;
clc;	clear;	close all;	
whitebg('k');	colormap iris;
plottools('off')

% Supress the warning raised by imshow() because why does this even exist?
warning('off','images:imshow:magnificationMustBeFitForDockedFigure')

%% Load image and begin processing

vFilenames = { 'Roman_Font_Test_Small.tif' };
vFilenames = { 'Roman_Font_Test.tif' };
% vFilenames = { 'Rosetta_Stone_All_Test.tif' };
% vFilenames = { 'Rosetta_Stone_All.tif' };
%vFilenames = { 'Rosetta_Stone_All.tif' 'Roman_Font_Test.tif' };

disp(sprintf('\n\nRestarting....\n\n'))
mGlyphPolygonsAll = [];

for k = 1:length(vFilenames)

	strFilename = vFilenames{k};
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
	%figure(1);
	%imshowz(label2rgb(imgGlyphs, @iris, [.5 .5 .5]));
	%hold on;

	nPts = 100;

	mGlyphPolygons = zeros(nPts*2,nGlyphs);
	for i = 1:length(vGlyphs)
		vGlyphPolygon = vGlyphs{i};
		%plot(iWidth-vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 2);

		% Interpolate the polygon so that all glyphs have the same # of points
		vGlyphPolygon = interppoly(vGlyphPolygon, nPts);

		% Because the image was flipped, the polygon x values will be wrong
		vGlyphPolygon(:,1) = iWidth-vGlyphPolygon(:,1);

		% Use this plot function when the image is transposed
		%plot(vGlyphPolygon(:,1), vGlyphPolygon(:,2), 'w', 'LineWidth', 1);
		% Use this plot to leave the original image unchanged
		%plot(vGlyphPolygon(:,2), vGlyphPolygon(:,1), 'w', 'LineWidth', 1);

		% Subtract the centroid to zero out polygons and allow analysis
		vGPCenter = mean(vGlyphPolygon);
% 		for j = 1:2
% 			vGlyphPolygon(:,j) = vGlyphPolygon(:,j)-vGPCenter(j);
% 		end
		vGlyphPolygon = vGlyphPolygon - ones(nPts,1)*vGPCenter;
		
		% Scale each polygon to have an average distance from center of 1
		fMeanDist = mean( sqrt(sum( vGlyphPolygon.^2,2)) );
		vGlyphPolygon = vGlyphPolygon/fMeanDist;

		%text(vGPCenter(1),vGPCenter(2),num2str(i),'Color','k','BackgroundColor','w');

		% Add the glyphs polygons to a big matrix
		mGlyphPolygons(:,i) = vGlyphPolygon(:);
	end

	% Transpose mGlyphPolygons so that they're oriented properly 
	% [Glyphs Points]
	mGlyphPolygons = mGlyphPolygons';
	
	mGlyphPolygonsAll = [ mGlyphPolygonsAll ; mGlyphPolygons ];
end

mGlyphPolygons = mGlyphPolygonsAll;

%%



%% Get Principal Components and plot a subset of axes
[mTPC,mGPPC,~] = getPC(mGlyphPolygons);
%mTPC = eye(size(mTPC,1));
mGPPC = mGlyphPolygons * mTPC;


figure(31);
clf;
hold on;
	
nSteps = 20;	
vColors = iris(nSteps);
for iIter = 1:nSteps
	
	

	%%% Use some subset of PCs to create a clustering tree

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

	vDistMetrics = {
		'sqeuclidean'
		'cityblock'
		'cosine'
		'correlation'
		'hamming' 
	};

	j = 3;

	nDims = size(mGPPC, 2);
	mGPPC = mGPPC(:,1:nDims);

	if( size(mGPPC,1) < 2 )
		
	end

	%%%%% Cluster

	if( size(mGPPC,1) < 2 )
		break;
	end
	
	nClusters = min( 2, size(mGPPC,1) )

	[T, mCentroids] = kmeans( mGPPC, nClusters, 'Distance', vDistMetrics{j} );
	%clc



	%%% Plot clusters

	% Sort clusters according to number of members and display the biggest
	% first
	nT = histc(T,1:nClusters);
	[nT,iSortedClusters] = sort(nT,'descend');

	iLargestCluster = iSortedClusters(1);

% 	mGlyphPolygons = mGlyphPolygons(T==iLargestCluster,:);
	mGPPC = mGPPC(T==iLargestCluster,:);
	
	x = mGPPC(:,1);
	y = mGPPC(:,2);
	z = mGPPC(:,3)*0 + iIter; 
	
	mGPCluster = mGlyphPolygons(T==iLargestCluster,:);
	vPolygon = {};
	nGPCluster = size(mGPCluster,1);
	
	fScale = 0.5;%iIter;
	 
	for j = 1:nGPCluster
		vPolygon{j} = reshape( mGPCluster(j,:), [], 2 );
		vPolygon{j} = [ vPolygon{j} vPolygon{j}(:,1)*0 ];
		vPolygon{j} = vPolygon{j} * fScale;
		
		
		fill3(vPolygon{j}(:,1)+x(j), ...
			vPolygon{j}(:,2)+y(j), ...
			vPolygon{j}(:,3)+z(j), vColors(iIter,:), ...
									'EdgeAlpha', 0.8);
		
		set(gca, 'ydir', 'reverse');
		rotate3d on;
	end
	
	 

	drawnow;
	
end

return;

iStart = find(nT < mean(nT), 1, 'first');
iEnd = find( nT < median(nT), 1, 'first' )+1;
iStart = 1;
iEnd = nClusters;d
for i = iSortedClusters(iStart:iEnd)'
	n = sum(T==i)
	if( n > 1 )
		figure(i);
		ShowGlyph(mGlyphPolygons, find(T==i));
		title(i);
		getframe;
	end
	%[ i n ]
end

sprintf('%1.4f - %d: %s', skewness(T), j, vDistMetrics{j} )
%pause




return;




%% K Nearest Neighbor Classifier
vDistMetrics = {
	'euclidean'
	'cityblock'
	'cosine'
	'correlation'
	'hamming' 
};

iNN = 10;	%ceil(sqrt(mean(nT)));
Mdl = fitcknn( mGPPC, T,'NumNeighbors', iNN, 'Distance', vDistMetrics{j} );

label = predict(Mdl,mGPPC);


plot(T+randn(size(T))/10,label+randn(size(T))/10,'.')
sum(T==label)/length(T)



%%% Get input

% close all;
% figure(1);
% plot([-1 101],[-1 101],'.');
% axis([0 100 0 100]);
% axis off;
% P = get_pencil_curve;

img = frame2im( getframe(gcf) );
img = img(:,:,1);
img = (img > 128);
img = imdilate( ~img, strel( 'disk', 20 ) ); 

figure(11);
imshow(img);

%%%


[vGlyphs,imgGlyphs] = bwboundaries((img)','noholes');
vGlyphPolygon = vGlyphs{1};
%vGlyphPolygon(:,1) = max(vGlyphPolygon(:,1))-vGlyphPolygon(:,1);
%vGlyphPolygon(:,2) = max(vGlyphPolygon(:,2))-vGlyphPolygon(:,2);
vGlyphPolygon = interppoly(vGlyphPolygon, nPts);

% Translate each polygon so that it is centered at the origin
vGPCenter = mean(vGlyphPolygon);
vGlyphPolygon = vGlyphPolygon - ones(nPts,1)*vGPCenter;

% Scale each polygon to have an average distance from center of 1
fMeanDist = mean( sqrt(sum(vGlyphPolygon.^2,2)) );
vGlyphPolygon = vGlyphPolygon/fMeanDist;

figure(21);
plotxy(vGlyphPolygon);
hold on;
plotxy(vGlyphPolygon(1,:),'or');
hold off;
set(gca,'ydir','reverse');
		
		
%%%

mGPPCTest = vGlyphPolygon(:)' * mTPC;
mGPPCTest = mGPPCTest(:,1:nDims);
label = predict(Mdl,mGPPCTest)
ShowGlyph(mGlyphPolygons, find(T==label));

% Take the distances between PC representations of glyphs and test glyph
vDist = sum( (ones(size(mGPPC,1),1)*mGPPCTest - mGPPC).^2, 2 );
% Sort them descending so that the smallest distaces are at the bottom
% This makes them easily visible in the command window
% [vDist, vSortIndex] = sort( vDist, 'descend' );
% vLabels = T(vSortIndex)
% vLabelsUnique = unique( vLabels(max(1,end-10):end) );
% 
% for i = 1:length(vLabelsUnique)
% 	figure(vLabelsUnique(i)+200);
% 	ShowGlyph(mGlyphPolygons, find(T==vLabelsUnique(i)));
% 	title(vLabelsUnique(i));
% end


