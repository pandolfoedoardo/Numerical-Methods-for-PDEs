function [Aloc,Bloc,Cloc,Area]=localBasis(Nelem,triang,coord)
% This function computes the local gradients of the linear shape functions
% (dN/dx and dN/dy) and the signed area for each triangular element.
% -- INPUTS --   
%        Nelem : number of elements (rows in triang)
%        triang: Nelem x 3 array of node indices (each row lists the 3 vertices of a triangle)
%        coord : Nnodes x 2 array of node coordinates; coord(i,:) = [x_i, y_i]
% -- OUTPUTS --
%       Aloc,Bloc,Cloc   : 3 x Nelem matrix with the linear coefficient of the planes we want to find
Aloc = zeros(3, Nelem);
Bloc = zeros(3, Nelem);
Cloc = zeros(3, Nelem);
Area = zeros(1, Nelem);
% Build the Vandermonde matrix
V=ones(3,3);
% Loop over elements
for i = 1:Nelem
    % element node indices, vertices of the triangles
    idx = triang(i,:);
    % Coordinates (x,y) of the vertices of the triangle 
    x = coord(idx,1); 
    y = coord(idx,2);
    % Build the matrix V for the element i V=[1,x,y]
    V(:,2)=x;
    V(:,3)=y;
    % Compute the area
    Area(i) =det(V)/2;
    % compute the coefficients Aloc Bloc, Cloc by solving the system V C= id
    % the plane we are looking for will have equation z= C(1,j)+C(2,j)x+C(3,j)y
    C=V\eye(3);
    Aloc(:,i)=C(1,:);
    Bloc(:,i)=C(2,:);
    Cloc(:,i)=C(3,:);
end
end