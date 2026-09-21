function [SysMat, rhs] = imposeBCDir(NDir, DirNod, DirVal, SysMat, rhs)

    
    Nodes = length(rhs);

    % STEP 1: Creazione del vettore di Lifting (u_g)
    % Creiamo un vettore lungo quanto tutti i nodi, pieno di zeri.
    u_g = zeros(Nodes, 1);
    
    % Inseriamo i valori noti SOLO nei nodi di bordo (i nostri DirVal)
    u_g(DirNod) = DirVal;
    
 %modfica del termine noto a destra
    rhs = rhs - SysMat * u_g; %è la parte destra tutto in una 
    
%imponiamo 1*Ui=gi modificando la SysMAt
    
    for i = 1:NDir
        n = DirNod(i); % Indice del nodo di bordo corrente
        
        % 3A. Azzera l'intera riga 'n' e colonna 'n'
        SysMat(n, :) = 0;
        SysMat(:, n) = 0;
        
        % 3B. Metti un bel 1 sulla diagonale principale
        SysMat(n, n) = 1;
        
        % 3C. Forza il termine noto ad assumere esattamente il valore esatto
        rhs(n) = DirVal(i);
    end
    %in questo modo abbiamo proprio 1 * U_bordo = Valore_bordo.
end