clear all; close all; clc;
% Define the 2D mesh and equation BC's and coefficients
meshdir = 'mesh-conv-diff'; % Imposta qui la singola mesh che vuoi analizzare
[Nelem, Nodes,triang, coord, NDir, DirNod,DirVal,force]=input_data(meshdir);
beta = [1,3];
Diff=1e-3;
% calculate element area and elemental coefficients of basis functions (b,c)
[Aloc,Bloc,Cloc,Area]=localBasis(Nelem,triang,coord);
%build stiffness matrix (without BCs)
stiffMat=stiffBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,Diff);
[rhs] = forcing2d(Nelem, Nodes, triang, coord, Area, force);
% build convection matrix (without BCs)
convMat=convBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,beta);
sysMat1 = stiffMat + convMat;
rhs1= rhs;
[sysMat1,rhs1]=imposeBCDir(NDir,DirNod,DirVal,sysMat1,rhs1);
uh1=sysMat1\rhs1;
%----------PLOT RISULTATO------------%
figure('units', 'normalized', 'outerposition', [0 0 1 1], 'Name', 'Galerkin vs SUPG');
subplot(1,2,1)
trisurf(triang,coord(:,1),coord(:,2),uh1);
title('Galerkin standard (no delta)');
view(45, 30);
    
%build stabilization matrix (without BCs)
delta = 5e-2; %ATTENZIONE se rifinissimo la mesh dovrei scrivere delta = Tau*Peclet number per avere consistenza
stabMat=stabBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,delta,beta);
    
sysMat=stiffMat + convMat +stabMat;

%----------DIRICHLET------------%

[sysMat,rhs]=imposeBCDir(NDir,DirNod,DirVal,sysMat,rhs);
uh=sysMat\rhs;
%----------PLOT RISULTATO------------%
subplot(1,2,2)
trisurf(triang,coord(:,1),coord(:,2),uh);
title(sprintf(' SUPG stabilized (\\delta = %g)', delta));
view(45, 30);
% plot_tri(coord,triang,1,'node_index_l')
% plot_tri(coord,triang,1,'edges')
% plot_tri(coord,triang,1,'face_index')