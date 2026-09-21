function [Nelem,Nodes,triang,coord,NDir,DirNod,ureal,force,u_0]=...
    input_data(meshdir)

%so che devo gestire in base a com'è la mesh
meshtriang=strcat(meshdir, '/triang.dat');
meshcoord=strcat(meshdir, '/xy.dat');
meshDirNod=strcat(meshdir, '/dirnod.dat'); %for Dirichlet Nodes
meshNeuNod = strcat(meshdir, '/neunod.dat'); %for Neumann Nodes
%meshDirVal=strcat(meshdir, '/dirval.dat'); 

load(meshtriang);
load(meshcoord);
load(meshDirNod);
load(meshNeuNod);
%load(meshDirVal); 

triang=triang(:,1:3);
coord=xy(:,1:2);
DirNod=dirnod; %change it dependinig on the nodes
NeuNod=neunod;
%DirVal=dirval;

DirNod = unique([dirnod; neunod]); 

% Ora non ci sono più nodi di Neumann in tutto il problema!
NeuNod = [];

Nelem=size(triang,1);
Nodes=size(coord,1);
NDir=length(DirNod);

ureal=@(x,y,t) exp(-t).*sin(pi.*x).*cos(pi.*y/2);
force=@(x, y, t) ((5 * pi^2 / 4) - 1) * exp(-t) .* sin(pi * x) .* cos((pi / 2) * y);
u_0=@(x,y) sin(pi.*x).*cos(pi.*y/2);
% Definizione della forzante f(x, y, t)

end
