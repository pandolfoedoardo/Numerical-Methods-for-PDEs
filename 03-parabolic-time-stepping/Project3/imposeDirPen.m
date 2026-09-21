function [stiffMat,rhs]=imposeDirPen(NDir,DirNod,DirVal,stiffMat,rhs)
%this function imposes the boundary conditions, using a "cheat": the
%penalty method we are considering Au=b where b is rhs and A the stiffmat
%here we are saying 1e15*u1+ A_1,2*u2+...+A_1,nodes*u_Nodes=1e15*b1 which
%approximatively implies u1=b1 hence giving the value we want to the
%boundary. And we do this for all the DirNodes
penalty= 1e15;
for i = 1:NDir
    stiffMat(DirNod(i),DirNod(i))= penalty;
    rhs(DirNod(i))=penalty*DirVal(i); 
end