function rhs=imposeBCNeu(NNeu_edges, NeuEdges, g_Neu, coord, rhs)
    % IMPOSEBCNEU impone le condizioni al contorno di Neumann.
    % le matrici di rigidezza e
    % di massa non vengono alterate. Si agisce solo sul termine noto.
    
    % Cicliamo sui lati del bordo che ci interessano
    for i = 1:NNeu_edges
        % Estraiamo i due nodi del segmento
        n1 = NeuEdges(i, 1);
        n2 = NeuEdges(i, 2);
        
        % Coordinate
        x1 = coord(n1, 1); y1 = coord(n1, 2); %in realtà x è sempre 1 o -1, non servirebbe
        x2 = coord(n2, 1); y2 = coord(n2, 2);
        
        % Lunghezza esatta del segmento (Pitagora)
        L = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        % Valore esatto del flusso sui due nodi
        g1 = g_Neu(x1, y1);
        g2 = g_Neu(x2, y2);
        
        % Regola dei Trapezi 1D esatta: g * L / 2, contando che usiamo base
        % di lagrange
        rhs(n1) = rhs(n1) + g1 * L / 2;
        rhs(n2) = rhs(n2) + g2 * L / 2;
    end
end
    

