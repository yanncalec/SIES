function IM=loadimgfnc(FileName)
% IM=loadimgfnc(FileName) 
%
% This function is borrowed from the code of Y.Capdeboscq and
% A.B.Karrman. It loads and binarizes an image. The output is a
% sparse binary matrix indicating the position of the boundary.

% We turn the image into a matrix of binary values:
IM = imread(FileName);

if length(size(IM)) == 3
    IM = flipud(im2double(rgb2gray(IM))) ;
else
    IM = flipud(im2double(IM)) ;
end
IM = round((IM-max(max(IM)))/min(min(IM-max(max(IM))))) ;

if max(size(IM)) > 500
    error('Size error: maximum size of the image is 500X500');
else
    % We use this reference value to find which values include the shape
    % and which do not:
    [rowPART,colPART] = find(IM == 1) ;
    [rowNOTPART,colNOTPART] = find(IM == 0) ;

    % We make the matrix holding the image sparse and then trim the edges
    % and add a 3 pixels to the edges
    IM = sparse(rowNOTPART,colNOTPART,1) ;
    IM = IM(min(rowPART):max(rowPART),min(colPART):max(colPART)) ;
    IMtrimmed = ones(size(IM)+6) ;
    IMtrimmed(4:end-3,4:end-3) = full(IM) ;
    IM = sparse(IMtrimmed) ;

    % We plot a contour plot of the image matrix:        
    % contour(IM) ;
end  
