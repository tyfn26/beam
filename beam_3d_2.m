clearAllMemoizedCaches;
clear all;
clc;
clearvars;

koordinat=[0 0 0; %1
          100 0 0;
          100 100 100]; %2
for jj=1:size(koordinat,1)-1
  conn_ilk(jj,1)=jj;
  conn_ilk(jj,2)=jj+1;
endfor

%conn_ilk=[1 2 ; 2 3 ];
eleman_sayisi=size(conn_ilk,1);
radius_number=eleman_sayisi-1;
radius_mesh=200;
rad=19.05;
N=200;
rho=7850e-9; %kg/mm3
node=zeros(0,3);
if eleman_sayisi==1
  node=linspace(koordinat(conn_ilk(1,1),:),koordinat(conn_ilk(1,2),:),N)';
else
    for j=1:radius_number
      [T1(j,:),T2(j,:),R]=radius(koordinat(conn_ilk(j,1),:),koordinat(conn_ilk(j,2),:),koordinat(conn_ilk(j+1,1),:),koordinat(conn_ilk(j+1,2),:),rad,radius_mesh);
      if radius_number<2
        node=vertcat(node,linspace(koordinat(conn_ilk(j,1),:),T1(j,:),N)');
        node=vertcat(node,R);
        node=vertcat(node,linspace(T2(j,:),koordinat(conn_ilk(j+1,2),:),N)');
      else
        if j<2
          node=vertcat(node,linspace(koordinat(conn_ilk(j,1),:),T1(j,:),N)');
          node=vertcat(node,R);
        else
          node=vertcat(node,linspace(T2(j-1,:),T1(j,:),N)');
          node=vertcat(node,R);
        endif
        if j==radius_number
          node=vertcat(node,linspace(T2(j,:),koordinat(conn_ilk(j+1,2),:),N)');
        endif
      endif
  endfor
endif
fprintf('Mesh atıldı. Nodelar birleştiriliyor...\n');
for k=1:size(node,1)-1
  conn(k,1)=k;
  conn(k,2)=k+1;
endfor
fprintf('Nodelar birleştirildi. Global stiffness matrisi oluşturuluyor...\n');
D=12.7; t=0.7;
di=D-2*t;
poisson=0.28;
%kappa=(4+3*poisson)/(2*(1+poisson));
kappa=2*(1+poisson)/(4+poisson); %shear correction factor(https://doi.org/10.18280/ti-ijes.630105)
%kappa=2; %shear correction factor
A=3.14*((D^2)-(di^2))/4;
E=204999.9984; %MPa
G=79999.99987; %MPa
I=3.14*((D^4)-(di^4))/64;
J=3.14*((D^4)-(di^4))/32;
nn=size(node,1); %number of nodes
ne=size(conn,1); %number of elements
ndof=nn*6;
K=sparse(ndof,ndof);
M=sparse(ndof,ndof);
f=zeros(ndof,1);
d=zeros(ndof,1);
Sigma=zeros(6*ne,1);

BC1=6*nn-5; %X
BC2=6*nn-4; %Y
BC3=6*nn-3;
BC4=6*nn-2;
BC5=6*nn-1;
BC6=6*nn;

%f(BC1)=100;
f(BC2)=100;
%f(BC3)=-100;
%f(BC4)=10000;
%f(BC5)=10000;
%f(BC6)=10000;

ifix=(1:6);
%ifix=horzcat(ifix,501:506);
%ifix=horzcat(ifix,1701:1706);

tip_dof=BC1:BC6;
isolve=setdiff(1:ndof,ifix);
free_dof=setdiff(1:ndof,[ifix,BC1:BC6]);

for e=1:ne
    n1=conn(e,1);  % element connectivity
    n2=conn(e,2);
    sctr=[ 6*n1-5 6*n1-4 6*n1-3 6*n1-2 6*n1-1 6*n1 6*n2-5 6*n2-4 6*n2-3 6*n2-2 6*n2-1 6*n2 ];
    K(sctr,sctr)=K(sctr,sctr) + kmat( node(conn(e,:),:),E,G,I,J,A,kappa);
    M(sctr,sctr)=M(sctr,sctr) + mmat( node(conn(e,:),:),rho,A,J,I);
end
fprintf('Global stiffness matrisi oluşturuldu.Matris çözümü gerçekleştiriliyor...\n');
if max(abs(d(BC1:BC6)))>0;
  K_ff=K(free_dof,free_dof);
  K_fp=K(free_dof,tip_dof);
  K_pf=K(tip_dof,free_dof);
  K_pp=K(tip_dof,tip_dof);
  d(free_dof)=K_ff \ (-K_fp * d(tip_dof));
  f(tip_dof)=K_pf*d(free_dof) + K_pp * d(tip_dof);
  %f(ifix)=K(ifix,:)*d;
  f=K*d;
else
  d(isolve)=K(isolve,isolve)\f(isolve);
  f=K*d;
endif
fprintf('Matris çözümü gerçekleştirildi. Reaksiyon kuvvetleri hesaplandı.\n');
for e=1:ne
  coord=node(conn(e,:),:);
  n1=coord(1,1:3); n2=coord(2,1:3);
  [T,TT]=transformation(coord);
  L=sqrt(dot((n2-n1),(n2-n1)));
  fi=(12*E*I)/(kappa*G*A*L^2);
  fi_=1/(1+fi);
  u1=6*e-5; %1
  v1=6*e-4; %2
  w1=6*e-3; %3
  tx1=6*e-2;%4
  ty1=6*e-1;%5
  tz1=6*e;%6
  u2=6*e+1;%7
  v2=6*e+2;%8
  w2=6*e+3;%9
  tx2=6*e+4;%10
  ty2=6*e+5;%11
  tz2=6*e+6;%12
  d_L=[ d(6*e-5) d(6*e-4) d(6*e-3) d(6*e-2) d(6*e-1) d(6*e) d(6*e+1) d(6*e+2) d(6*e+3) d(6*e+4) d(6*e+5) d(6*e+6)]';
  d_L=TT*d_L;
  Sigma(6*e-4)=E*D*(d_L(12)-d_L(6))/(2*L);
  %Sigma(6*e-4)=E*t*(d_L(12)-d_L(6))/L;
  Sigma(6*e-3)=E*D*(d_L(11)-d_L(5))/(2*L);
  Sigma(6*e-5)=(E*(d_L(7)-d_L(1))/L)+sqrt(Sigma(6*e-4)^2 + Sigma(6*e-3)^2);
  %Sigma(6*e-5)=(E*(d_L(7)-d_L(1))/L);
  Sigma(6*e-2)=G*D*(d_L(10)-d_L(4))/(2*L);
  Sigma(6*e-1)=-G*fi*fi_*(2*d_L(2)+d_L(6)*L-2*d_L(8)+d_L(12)*L)/(2*L);
  Sigma(6*e)=-G*fi*fi_*(2*d_L(3)-d_L(5)*L - 2*d_L(9)-d_L(11)*L)/(2*L);


  tensor=[Sigma(6*e-5) Sigma(6*e-2) Sigma(6*e);
          Sigma(6*e-2) Sigma(6*e-4) Sigma(6*e-1);
           Sigma(6*e) Sigma(6*e-1) Sigma(6*e-3)];
  [,mp]=eig(tensor);
  max_p(e,:)=max(mp);
end
fprintf('Stress değerleri hesaplandı.\n');
for n=1:nn
  dx(n)=d(6*n-5);
  dy(n)=d(6*n-4);
  dz(n)=d(6*n-3);
  dtx(n)=d(6*n-2);
  dty(n)=d(6*n-1);
  dtz(n)=d(6*n);
endfor
fprintf('Nodelara deplasmanlar işlendi.\n');
dx_max=max(dx);
dy_max=max(dy);
dz_max=max(dz);
dtx_max=max(dtx);
dty_max=max(dty);
dtz_max=max(dtz);
deformed_X=node(:,1)+dx';
deformed_Y=node(:,2)+dy';
deformed_Z=node(:,3)+dz';

##figure;
##hold on;
##plot3(node(:,1),node(:,2),node(:,3),'LineWidth', 1);
##%plot3(deformed_X,deformed_Y,deformed_Z,'r','LineWidth', 1);
##view(3);
##X=node(size(node,1),1);
##Y=node(size(node,1),2);
##Z=node(size(node,1),3);
##U=f(BC1); V=f(BC2); W=f(BC3);
##quiver3(X,Y,Z,U,V,W,'r');
##
##axis equal;
##grid on;
##xlabel('X');
##ylabel('Y');
##zlabel('Z');
##legend('Undeformed', 'Deformed');
%clc;
fprintf('Analiz tamamlandı.\n');
disp_x=dx(n)
disp_y=dy(n)
disp_z=dz(n)
max_principal=max(max_p)

    %% 5. Solve the Eigenvalue Problem
    % Solve: K * Phi = Lambda * M * Phi
    numModes = 12; % Number of modes to compute
    [Phi, Lambda] = eigs(K, M, numModes, 'sm');

    % Extract eigenvalues and convert to natural frequencies
    lambda = diag(Lambda);
    omega_n = sqrt(lambda);   % Natural frequencies in rad/s
    freq_Hz = omega_n / (2*pi); % Natural frequencies in Hz

    disp('Natural Frequencies (Hz):');
    disp(freq_Hz);

    %% 6. Post-Processing: Reconstruct full mode shapes
    % Initialize a matrix for full mode shapes (including fixed DOF)
    fullModeShapes = zeros(totalDof, numModes);
    for i = 1:numModes
        % Place the computed eigenvectors into the full vector
        fullModeShapes(freeDof, i) = Phi(:, i);
        % The fixed DOF remain zero (default)
    end

    %% 7. Visualize the Mode Shapes
    % Plot the undeformed beam
    figure('Position', [100, 100, 1200, 800], 'Name', '3D Timoshenko Beam Mode Shapes');
    subplot(2, 3, 1);
    plot3(nodeCoords, zeros(numNodes,1), zeros(numNodes,1), 'k-o', 'LineWidth', 2, 'MarkerFaceColor', 'k');
    title('Undeformed Beam');
    xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
    axis equal; grid on;

    % Plot each mode shape
    scaleFactor = 0.2; % Scaling factor for deformation visualization
    for mode = 1:numModes
        subplot(2, 3, mode+1);
        hold on;

        % Extract displacements for this mode
        U = fullModeShapes(1:dofPerNode:end, mode);        % u_x
        V = fullModeShapes(2:dofPerNode:end, mode);        % u_y
        W = fullModeShapes(3:dofPerNode:end, mode);        % u_z

        % Plot the deformed shape
        deformedX = nodeCoords + scaleFactor * U;
        deformedY = 0 + scaleFactor * V; % Add deformation to initial Y=0 position
        deformedZ = 0 + scaleFactor * W; % Add deformation to initial Z=0 position

        plot3(deformedX, deformedY, deformedZ, 'r-o', 'LineWidth', 2, 'MarkerFaceColor', 'r');
        plot3(nodeCoords, zeros(numNodes,1), zeros(numNodes,1), 'k--'); % Undeformed

        title(sprintf('Mode %d: %.2f Hz', mode, freq_Hz(mode)));
        xlabel('X (m)'); ylabel('Y (m)'); zlabel('Z (m)');
        axis equal; grid on;
        view(30, 30); % Set a nice 3D view
        hold off;
    end
