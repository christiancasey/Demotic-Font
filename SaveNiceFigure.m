function [ img ] = SaveNiceFigure( fScale, strFolder )
% SAVENICEFIGURE  Save a nice-looking version of the current figure.
% 
	if( exist('fScale') ~= 1 )
		fScale = 2.5;
	end
	
	if( exist('strFolder') ~= 1 )
		strFolder = 'Images';
	end
	
	[ vFiles, nFiles ] = getFileList( [strFolder '/*.png'] );
	
	strFile = sprintf('%s/%04i.png', strFolder, nFiles+1);
	strEval = sprintf('export_fig %s -m%2.1f', strFile, fScale);
	
	eval(strEval);
	
	img = imread(strFile);
end