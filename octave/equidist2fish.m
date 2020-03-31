pkg load image

equi = im2double(flipud(imread('Downloads/earth.jpg')));
printf('Image loaded.\n');
fflush(stdout);
N=256;
equi = imresize(equi, [N,N]);
%figure; imshow(equi);
printf('Image resized.\n');
fflush(stdout);

% for remap style mapping, must do the inverse transform,
% from fisheye(X,Y) to equi(x,y)

aperture = pi;

anglexrad = 0;

% normalizing to [-1, 1]
xfishr = 2.*(1:N) / N .- 1;
yfishc = (2.*(1:N) / N .- 1)';
for i=1:N
 xfish(i,:) = xfishr;
 yfish(:,i) = yfishc;
end
rfish = sqrt(xfish.*xfish + yfish.*yfish);
theta = atan2(yfish, xfish) + anglexrad;
phi = rfish*aperture/2;

% here, phi = 0 @ Pz = 1
Px = sin(phi).*cos(theta);
Py = sin(phi).*sin(theta);
Pz = cos(phi);

angleyrad = 70*pi/180;
angle2rad = 90*pi/180;


%rotationmatrixy = [ cos(-angle2rad), 0, sin(-angle2rad)
% 0, 1, 0;
% -sin(-angle2rad), 0, cos(-angle2rad) ]; 
% 
%rotationmatrixz = [cos(angle2rad), -sin(angle2rad), 0;
% sin(angle2rad), cos(angle2rad), 0;
% 0, 0, 1];
% 
% rotationmatrixx = [ 1, 0, 0;
% 0, cos(angleyrad), -sin(angleyrad);
% 0, sin(angleyrad), cos(angleyrad)]; 

% for i=1:N
% for j = 1:N
%  
%  Inputm = [Px(i,j); Py(i,j); Pz(i,j)];
%  Outputm = rotationmatrixx * Inputm;
%  Px(i,j) = Outputm(1,1); 
%  Py(i,j) = Outputm(2,1); 
%  Pz(i,j) = Outputm(3,1);
% 
% end
%end 
  PyR = cos(angleyrad) .* Py - sin(angleyrad) .* Pz;
	PzR = sin(angleyrad) .* Py + cos(angleyrad) .* Pz;
  PxR = Px;
  
  printf('Rotated around x. \n');
  fflush(stdout);
  
  PxR2 = cos(angle2rad).* PxR - sin(angle2rad) .* PyR;
  PyR2 = sin(angle2rad).* PxR + cos(angle2rad) .* PyR;
  PzR2 = PzR;
  
  printf('Rotated around z. \n');
  fflush(stdout);


longi 	= atan2(PyR2, PxR2);
lat	 	  = atan2(PzR2, sqrt(PxR2.*PxR2 + PyR2.*PyR2));
X = longi / pi;
Y = 2*lat / pi;
%X = Px;% this gives back a slightly distorted equi turned on its side
%Y = Py;
% these are [-1,1] normalized co-ords
Xindex = floor(X*N/2 + N/2 + 1);
Yindex = floor(Y*N/2 + N/2 + 1);

printf('Found indices. \n');
fflush(stdout);

%figure; surf(X,Y,  "CData", equi, "FaceColor", "texturemap");
%figure; surf(X,  "CData", equi, "FaceColor", "texturemap");
%figure; surf(Y,  "CData", equi, "FaceColor", "texturemap");

fishe = zeros(size(equi));

for i=1:N
 for j=1:N
  if(Xindex(i,j)<=N && Yindex(i,j)<=N && Xindex(i,j)>0 && Yindex(i,j)>0)
   fishe(i,j,:)=equi(N-Yindex(i,j)+1,Xindex(i,j),:);
  end
 end
end

imshow(fishe)