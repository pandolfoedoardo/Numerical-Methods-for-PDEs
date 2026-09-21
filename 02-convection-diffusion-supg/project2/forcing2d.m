function rhs = forcing2d(Nelem, Nodes, triang, coord, Area, force)
    % FORCING Assembla il vettore dei termini noti (RHS) per il 2D
    % Utilizza il Metodo dei Trapezi (Quadratura sui vertici).
    
    % Inizializza il vettore colonna pieno di zeri
    rhs = zeros(Nodes, 1);
    
    for i = 1:Nelem
        % Estraiamo gli indici dei 3 nodi del triangolo corrente
        n1 = triang(i, 1); 
        n2 = triang(i, 2); 
        n3 = triang(i, 3);
        
        % Estraiamo le coordinate esatte (x,y) dei 3 vertici
        x1 = coord(n1, 1); y1 = coord(n1, 2);
        x2 = coord(n2, 1); y2 = coord(n2, 2);
        x3 = coord(n3, 1); y3 = coord(n3, 2);
        
        % METODO DEI TRAPEZI:
        % Valutiamo la funzione della forza esattamente sui 3 vertici
        f1 = force(x1, y1); 
        f2 = force(x2, y2); 
        f3 = force(x3, y3); 
        
        % Assembliamo i contributi nel vettore globale
        % (Valore della forza sul nodo * Area / 3)
        rhs(n1) = rhs(n1) + f1 * Area(i) / 3;
        rhs(n2) = rhs(n2) + f2 * Area(i) / 3;
        rhs(n3) = rhs(n3) + f3 * Area(i) / 3;
    end
end

