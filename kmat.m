function [ke] = kmat(coord, e, G, I, J, A, kappa)
%kappa=0.5;%shear correction factor

EI=e*I;  GJ=G*J; EA = e*A;
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
fi=12*e*I/(kappa*G*A*L^2);
fi_=1/(1+fi);
ke = [EA/L, 0, 0, 0, 0, 0, -(EA/L), 0, 0, 0, 0, 0;
    0, (12*EI*fi_)/L^3, 0, 0, 0, (6*EI*fi_)/L^2, 0, -((12*EI*fi_)/L^3),0,0, 0, (6*EI*fi_)/L^2;
    0, 0, (12*EI*fi_)/L^3, 0, -((6*EI*fi_)/L^2), 0, 0, 0, -((12*EI*fi_)/L^3), 0, -((6*EI*fi_)/L^2), 0;
    0, 0, 0, GJ/L, 0, 0, 0, 0, 0, -(GJ/L), 0, 0;
    0, 0, -((6*EI*fi_)/L^2), 0, ((4+fi)*EI*fi_)/L, 0, 0, 0, (6*EI*fi_)/L^2, 0,((2-fi)*EI*fi_)/L, 0;
    0, (6*EI*fi_)/L^2, 0, 0, 0, ((4+fi)*EI*fi_)/L, 0, -((6*EI*fi_)/L^2), 0, 0, 0, ((2-fi)*EI*fi_)/L;
    -(EA/L), 0, 0, 0, 0, 0, EA/L, 0, 0, 0, 0, 0;
    0, -((12*EI*fi_)/L^3), 0, 0, 0, -((6*EI*fi_)/L^2),0, (12*EI*fi_)/L^3, 0, 0, 0, -((6*EI*fi_)/L^2);
    0, 0, -((12*EI*fi_)/L^3), 0, (6*EI*fi_)/L^2, 0, 0, 0, (12*EI*fi_)/L^3, 0, (6*EI*fi_)/L^2, 0;
    0, 0, 0, -(GJ/L), 0, 0, 0, 0, 0, GJ/L, 0, 0;
    0, 0, -((6*EI*fi_)/L^2), 0, ((2-fi)*EI*fi_)/L, 0, 0,0, (6*EI*fi_)/L^2, 0, ((4+fi)*EI*fi_)/L, 0;
    0, (6*EI*fi_)/L^2, 0, 0, 0, ((2-fi)*EI*fi_)/L, 0, -((6*EI*fi_)/L^2), 0, 0, 0, ((4+fi)*EI*fi_)/L];
ke = T*ke*TT;
