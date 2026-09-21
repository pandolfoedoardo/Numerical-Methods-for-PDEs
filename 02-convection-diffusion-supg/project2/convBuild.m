function convMat = convBuild(Nelem, Nodes, triang, Bloc, Cloc, Area, beta)

    convMat = zeros(Nodes, Nodes);
    
    % Ciclo su tutti gli elementi (triangoli)
    for i = 1:Nelem
        
        % Ciclo sulle righe locali (funzione di Test, indice i)
        for iloc = 1:3
            iglob = triang(i, iloc);
            
            % Ciclo sulle colonne locali (funzione Trial, indice j)
            for jloc = 1:3
                jglob = triang(i, jloc);
                
                % Calcolo del termine convettivo locale
                % jloc governa i gradienti (Bloc, Cloc)
                % iloc governa l'integrale della funzione di test (Area/3)
                termine_convettivo = (beta(1)*Bloc(jloc, i) + beta(2)*Cloc(jloc, i)) * (Area(i) / 3);
                
                % Assemblaggio Globale (Attenzione a sommare!)
                convMat(iglob, jglob) = convMat(iglob, jglob) + termine_convettivo;
                
            end
        end
    end
end
