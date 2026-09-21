function [Nelem,Nodes,triang,coord,NDir,DirNod,DirVal,force]=...
    input_data(meshdir)

%so che devo gestire in base a com'è la mesh
meshtriang=strcat(meshdir, '/triang.dat');
meshcoord=strcat(meshdir, '/xy.dat');
meshDirNod=strcat(meshdir, '/dirnod.dat'); %for Dirichlet Nodes
meshDirVal=strcat(meshdir, '/dirval.dat'); 

load(meshtriang);
load(meshcoord);
load(meshDirNod);
load(meshDirVal); 

triang=triang(:,1:3);
coord=xy(:,1:2);
DirNod=dirnod; %change it dependinig on the nodes
DirVal=dirval;

Nelem=size(triang,1);
Nodes=size(coord,1);
NDir=length(DirNod);


force=@(x,y) 0;
end
