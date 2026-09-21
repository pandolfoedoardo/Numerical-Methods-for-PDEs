function [Nelem,Nodes,triang,coord,NDir,DirNod,NNeu,NeuNod,ureal,force]=...
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

Nelem=size(triang,1);
Nodes=size(coord,1);
NDir=length(DirNod);
NNeu=length(NeuNod);

ureal=@(x,y) sin(pi.*x).*cos(pi.*y/2);
force=@(x,y) 5*pi^2*ureal(x,y)/4; % avevamo messo 0 come mai? per verificare con una funzioen semplice
%1 a sinsitra e 0 a destra
end
