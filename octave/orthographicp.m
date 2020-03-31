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

PxW = sin(lat).*cos(longi);
PyW = sin(lat).*sin(longi);
PzW = cos(lat);

figure; surf(PxP, PyP, PzP, "CData", equi, "FaceColor", "texturemap");
view(2); axis off; axis square; title('Paul')

figure; surf(PxW, PyW, PzW, "CData", equi, "FaceColor", "texturemap");
view(2); axis off; axis square; title('Wolfram')

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
  PxProt(i,j) = Outputm(1,1); 
  PyProt(i,j) = Outputm(2,1); 
  PzProt(i,j) = Outputm(3,1);
  
 end
end

for i=1:N
 for j = 1:N
  
  Inputm = [PxW(i,j); PyW(i,j); PzW(i,j)];
  Outputm = rotationmatrixx * Inputm;
  PxWrot(i,j) = Outputm(1,1); 
  PyWrot(i,j) = Outputm(2,1); 
  PzWrot(i,j) = Outputm(3,1);
  
 end
end

figure; surf(PxProt, PyProt, PzProt, "CData", equi, "FaceColor", "texturemap");
view(2); axis off; axis square;title('Paul')

figure; surf(PxWrot, PyWrot, PzWrot, "CData", equi, "FaceColor", "texturemap");
view(2); axis off; axis square;title('Wolfram')


for i=1:N
 for j = 1:N
  
  Inputm = [PxProt(i,j); PyProt(i,j); PzProt(i,j)];
  Outputm = rotationmatrixy * Inputm;
  PxProt2(i,j) = Outputm(1,1); 
  PyProt2(i,j) = Outputm(2,1); 
  PzProt2(i,j) = Outputm(3,1);
  
 end
end

figure; surf(PxProt2, PyProt2, PzProt2, "CData", equi, "FaceColor", "texturemap");
view(2); axis off; axis square;
