
function stiffMat = stiffBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, Diff)
    stiffMat = zeros(Nodes, Nodes);
    
    for i = 1:Nelem
        for j = 1:3
            for k = 1:3
                % Indici globali usando triang(i, j) e triang(i, k)
                nj = triang(i, j);
                nk = triang(i, k);
                
                % Prodotto scalare dei gradienti (sto risolvendo
                % l'integrale)
                K_loc = Diff * (Bloc(j,i)*Bloc(k,i) + Cloc(j,i)*Cloc(k,i)) * Area(i);
                
                % Assemblaggio
                stiffMat(nj, nk) = stiffMat(nj, nk) + K_loc;
            end
        end
    end
end