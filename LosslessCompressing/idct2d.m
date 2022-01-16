function Y = idct2d(dctmat)
[Rows, Cols] = size(dctmat);
 for k = 1:Cols
     dcty(:,k) = idct1d(dctmat(:,k)')';
 end
 for k = 1:Rows
     Y(k,:) = idct1d(dcty(k,:));
 end