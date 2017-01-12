function ShowGlyph(mGlyphPolygons, g)
	
	%figure;
	clf;
	hold on;
	
	vColors = iris(length(g));
	for i = 1:length(g)
		
		vPolygon = mGlyphPolygons(g(i),:);
		vPolygon = reshape( vPolygon, [], 2 );
		vPolygon(:,1) = vPolygon(:,1) + (i-1)*max(mGlyphPolygons(:)*2);

% 		plotxy( vPolygon );
		fill( vPolygon(:,1), vPolygon(:,2), vColors(i,:), 'EdgeColor', vColors(i,:) );
% 		hold on;
% 		plotxy(vPolygon(1,:),'or');
% 		hold off;
		set(gca,'ydir','reverse');
	end
	
	hold off;
