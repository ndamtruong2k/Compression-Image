function Y = dct2d(mat)
[Rows, Cols] = size(mat);
 for k = 1:Rows
     dctx(k,:) = dct1d(mat(k,:));
 end
 for k = 1:Cols
     Y(:,k) = dct1d(dctx(:,k)')';
 end