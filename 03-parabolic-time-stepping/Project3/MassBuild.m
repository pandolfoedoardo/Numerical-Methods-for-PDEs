function MassMat = MassBuild(Nelem, Nodes,triang,Area)
MassMat=sparse(Nodes,Nodes);
for i=1:Nelem
    for j=1:3
        nj=triang(i,j);
        for k=1:3
            nk=triang(i,k);
            if nj==nk
                MassMat(nj,nk)= MassMat(nj,nk)+Area(i)/6;
            else
                MassMat(nj,nk)= MassMat(nj,nk)+Area(i)/12;
            end
              
        end
    end
end

            


