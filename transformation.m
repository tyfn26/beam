function [T,TT] = transformation(coord)
  Lambda=zeros(3,3);
x1=coord(1,1);
x2=coord(2,1);
y1=coord(1,2);
y2=coord(2,2);
z1=coord(1,3);
z2=coord(2,3);
L = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +(z2-z1)*(z2-z1));
if x1 == x2 && y1 == y2
  if z2 > z1
    Lambda = [0 0 1 ; 0 1 0 ; -1 0 0];
  else
    Lambda = [0 0 -1 ; 0 1 0 ; 1 0 0];
  end
  else
      CXx = (x2-x1)/L;
CYx = (y2-y1)/L;
CZx = (z2-z1)/L;
D = sqrt(CXx*CXx + CYx*CYx);
CXy = -CYx/D;
CYy = CXx/D;
CZy = 0;
CXz = -CXx*CZx/D;
CYz = -CYx*CZx/D;
CZz = D;
Lambda = [CXx CYx CZx ;CXy CYy CZy ;CXz CYz CZz];
endif


T = zeros(12);
T([1, 2, 3], [1, 2, 3]) = Lambda;
T([4,5,6], [4,5,6]) = Lambda;
T([7,8,9], [7,8,9]) = Lambda;
T([10,11,12], [10,11,12]) = Lambda;
TT = T';
