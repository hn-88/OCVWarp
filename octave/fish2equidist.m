pkg load image

equi = im2double(flipud(imread('Downloads/earthorthogr.png')));
disp('Image loaded.');
N=256;
equi = imresize(equi, [N,N]);
%figure; imshow(equi);
disp('Image resized.');

longir = (1:N)*pi/(N/2) .- pi;
latc = ((1:N)*pi/N .- pi/2)';
for i=1:N
 longi(i,:) = longir;
 lat(:,i)   =   latc;
end


PxP = cos(lat).*cos(longi);
PyP = cos(lat).*sin(longi);
PzP = sin(lat);
aperture = pi;
R = 2 .* atan2(sqrt(PxP.*PxP + PyP.*PyP), PzP) ./ aperture;
theta = atan2(PyP, PxP);
X =  R .* cos(theta);
Y =  R .* sin(theta);
% these are [-1,1] normalized co-ords
Xindex = floor(X*N/2 + N/2 + 1);
Yindex = floor(Y*N/2 + N/2 + 1);

%figure; surf(X,Y,  "CData", equi, "FaceColor", "texturemap");
%figure; surf(X,  "CData", equi, "FaceColor", "texturemap");
%figure; surf(Y,  "CData", equi, "FaceColor", "texturemap");

fishe = zeros(size(equi));

for i=1:N
 for j=1:N
  if(Xindex(i,j)<257 && Yindex(i,j)<257 && Xindex(i,j)>0 && Yindex(i,j)>0)
   fishe(i,j,:)=equi(Xindex(i,j),Yindex(i,j),:);
  end
 end
end

fishe = flipud(fishe);
fishe = fliplr(fishe);

imshow(fishe)