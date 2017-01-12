function mouseMove (object, eventdata)

% 	try

		global vPts;
		global bButtonDown;
		global img;
		[h w] = size(img);

		C = get (gca, 'CurrentPoint');
		x = C(1,1);
		y = C(1,2);

		if( bButtonDown )
			vPts = [ vPts ; x y ];
			x = max( min(x, 1), 0 );
			y = max( min(y, 1), 0 );
% 			x = max( min( round(x), 1), w );
% 			y = max( min( round(y), 1), h );
% 			img(y,x) = 1.0;
% 			imshow(img);

			plotxy(vPts);
		end
% 	catch
% 		bButtonDown = false;
% 	end
	%title(gca, ['(X,Y) = (', x, ', ',y, ')']);
	
	
	
	
	
	
	