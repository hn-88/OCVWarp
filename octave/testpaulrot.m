X=[-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ;
-1 -0.75 -0.5 -0.25 0 0.25 0.5 0.75 1.0 ];

Y=[1 1 1 1 1 1 1 1 1;
-0.75 -0.75 -0.75 -0.75 -0.75 -0.75 -0.75 -0.75 -0.75 ;
-0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 -0.5 ;
-0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 -0.25 ;
0 0 0 0 0 0 0 0 0 ;
0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 0.25 ;
0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 0.5 ;
0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 0.75 ;
1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 ];

FOV = pi;

longitude  = X.*pi;
latitude = Y.*pi/2;

Px = cos(latitude) .* cos(longitude);
Py = cos(latitude) .*sin(longitude);
Pz = sin(latitude);

angleyrad = 30*pi/180;

rotationmatrixy = [ cos(-angleyrad), 0, sin(-angleyrad)
 0, 1, 0;
 -sin(-angleyrad), 0, cos(-angleyrad) ]; 
rotationmatrixx = [ 1, 0, 0;
 0, cos(angleyrad), -sin(angleyrad);
 0, sin(angleyrad), cos(angleyrad)]; 
rotationmatrixz = [cos(angleyrad), -sin(angleyrad), 0;
 sin(angleyrad), cos(angleyrad), 0;
 0, 0, 1]; 
 
Pxrot = Px;
Pyrot = Py;
Pzrot = Pz;

for i=1:3
 for j = 1:3
  
  Inputm = [Px(i,j); Py(i,j); Pz(i,j)];
  Outputm = rotationmatrixy * Inputm;
  Pxrot(i,j) = Outputm(1,1); 
  Pyrot(i,j) = Outputm(2,1); 
  Pzrot(i,j) = Outputm(3,1);
  
 end
end
r = 2 .* atan2(sqrt(Px.*Px + Pz.*Pz ),Py ) ./ FOV;
theta = atan2(Pz , Px );

x = round(100.*r.* cos(theta) );
y = round(100.*r .* sin(theta) );

x = x./100;
y = y./100;

for i = 1:9
  for j = 1:9
    printf("(%1.2f, %1.2f) maps to ", X(j,i), Y(j,i))
    printf("(%1.2f, %1.2f) \n", x(j,i), y(j,i))
  end
  
end