function err=evalErrH1(Nelem, coord, triang, Area, Bloc, Cloc, uh, grad_ureal,err)
err= err^2;
err_semi = 0;
for i = 1:Nelem


% Estraiamo i nodi e le coordinate del triangolo
    n = triang(i,:);
    x = coord(n, 1);
    y = coord(n, 2);
       
        % 1. CALCOLO ERRORE SUI GRADIENTI (Seminorma H1)
        % Troviamo il baricentro del triangolo --> quadratura esatto con
        % costanti
    xb = (x(1) + x(2) + x(3)) / 3;
    yb = (y(1) + y(2) + y(3)) / 3;
        
        % Valutiamo il gradiente ESATTO nel baricentro
    grad_ex = grad_ureal(xb, yb); 
    % Calcoliamo il gradiente NUMERICO (costante sul triangolo)
        % la formula viene dal fatto che uh è comb lineare dei tre elementi
        % della base di lagrange
        du_dx_h = (uh(n(1))*Bloc(1,i) + uh(n(2))*Bloc(2,i) + uh(n(3))*Bloc(3,i)); 
        
        du_dy_h = (uh(n(1))*Cloc(1,i) + uh(n(2))*Cloc(2,i) + uh(n(3))*Cloc(3,i));
        
        % Aggiungiamo l'errore del gradiente (Area * differenza al quadrato)
    err_semi = err_semi + Area(i) * ((grad_ex(1) - du_dx_h)^2 + (grad_ex(2) - du_dy_h)^2);
end
err= sqrt(err+err_semi);
         
    
   

end