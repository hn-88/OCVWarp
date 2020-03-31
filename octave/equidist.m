pkg load image

equi = im2double(flipud(imread('Downloads/earth.jpg')));
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

angleyrad = -90*pi/180;
angle2rad = 90*pi/180;


rotationmatrixy = [ cos(-angle2rad), 0, sin(-angle2rad)
 0, 1, 0;
 -sin(-angle2rad), 0, cos(-angle2rad) ]; 
 
rotationmatrixz = [cos(angleyrad), -sin(angleyrad), 0;
 sin(angleyrad), cos(angleyrad), 0;
 0, 0, 1];
 
 rotationmatrixx = [ 1, 0, 0;
 0, cos(angleyrad), -sin(angleyrad);
 0, sin(angleyrad), cos(angleyrad)]; 

 for i=1:N
 for j = 1:N
  
  Inputm = [PxP(i,j); PyP(i,j); PzP(i,j)];
  Outputm = rotationmatrixx * Inputm;
  PxP(i,j) = Outputm(1,1); 
  PyP(i,j) = Outputm(2,1); 
  PzP(i,j) = Outputm(3,1);
  
 end
end




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
  if(Xindex(i,j)<=N && Yindex(i,j)<=N && Xindex(i,j)>0 && Yindex(i,j)>0)
   fishe(Xindex(i,j),Yindex(i,j),:)=equi(i,j,:);
  end
 end
end

imshow(fishe)