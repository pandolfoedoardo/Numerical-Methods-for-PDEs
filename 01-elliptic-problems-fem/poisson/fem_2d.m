%%% Finite element code for 2D piecewise linear Galekrin 
%%% -div(Diff grad u)(x)=f(x) with Dirichlet BC 
clear all; close all; clc;

% Define the 2D mesh and equation BC's and coefficients
Meshes = {'mesh0', 'mesh1', 'mesh2', 'mesh3', 'mesh4-renum'};
errors= zeros(length(Meshes),1);
h = zeros(length(Meshes),1);
for i= 1:5 
    meshdir=Meshes{i};
    [Nelem,Nodes,triang,coord,NDir,DirNod,NNeu, NeuNod, ureal,force]=input_data(meshdir);
    Diff=1;
    
    % calculate element area and elemental coefficients of basis functions
    % (b,c)
    [Aloc,Bloc,Cloc,Area]=localBasis(Nelem,triang,coord);
    
    %build stiffness matrix (without BCs)
    stiffMat=stiffBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,Diff);
    
    
    [rhs] = forcing2d(Nelem, Nodes, triang, coord, Area, force);
    
    % build convection matrix (without BCs)
    %convMat=convBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,beta);
    
    %build stabilization matrix (without BCs)
    %stabMat=stabBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,delta);
    
    % sysMat=stiffMat+convMat+stabMat;
    
    % impose BCs
    %[sysMat,rhs]=imposeBC(Nodes,NDir,DirNod,DirVal,sysMat);
        %----------NEUMANN-----------%

    %%% ATTENZIONE %%% nei nodi di dirichlet la condizione di dirichlet
    %%% spazza via quelal di neumann
    sysMat = stiffMat;
    g_Neu=@(x,y) -pi*sign(x).*cos(pi/2*y); %tiene in cosiderazione entrambi i lati
    NeuEdges = [];

    for k = 1:Nelem %è complicato perchè deve gestire i nodi agli angoli
        n = triang(k, :); % I tre nodi del triangolo [n1, n2, n3]
        
        % Definiamo le coppie di nodi che formano i tre lati
        lati = [n(1), n(2); 
                n(2), n(3); 
                n(3), n(1)];
                
        for j = 1:3
            nodeA = lati(j,1);
            nodeB = lati(j,2);
            
            % CONDIZIONE AGGIORNATA:
            % Il lato è Neumann se:
            % (A e B sono Neu) OPPURE (A è Neu e B è Dir) OPPURE (A è Dir e B è Neu)
            if (ismember(nodeA, NeuNod) && ismember(nodeB, NeuNod)) || ...
               (ismember(nodeA, NeuNod) && ismember(nodeB, DirNod)) || ...
               (ismember(nodeA, DirNod) && ismember(nodeB, NeuNod))
           
                % IMPORTANTE: In una mesh 2D, un lato di bordo appartiene a UN SOLO triangolo.
                % Questo garantisce che non aggiungeremo lati interni per errore, 
                % a patto che i nodi Dirichlet siano solo sul bordo.
                NeuEdges = [NeuEdges; nodeA, nodeB];
            end
        end
    end

    NNeu_edges = size(NeuEdges, 1);
    
    % 2. Chiamiamo la funzione corretta passando i lati e la function handle
    rhs = imposeBCNeu(NNeu_edges, NeuEdges, g_Neu, coord, rhs);

    
    %----------DIRICHLET------------%
    DirVal=ureal(coord(DirNod, 1), coord(DirNod, 2));
    [sysMat,rhs]=imposeDirPen(NDir,DirNod,DirVal,sysMat,rhs);%devo capire i Dirval

    uh=sysMat\rhs;
    errors(i)=evalErr(Nelem,coord,triang,Area,uh,ureal);
    h(i)=max(2*sqrt(Area));
    logh(i)=log(h(i));
    logerr(i)=log(errors(i));

    grad_ureal = @(x,y) [pi * cos(pi.*x) .* cos(pi.*y/2), -(pi/2) * sin(pi.*x) .* sin(pi.*y/2)];
    errorsH(i)=evalErrH1(Nelem, coord, triang, Area, Bloc, Cloc, uh, grad_ureal,errors(i)); %da capire bene se è migliorabile
    logerrH(i)=log(errorsH(i));

end
%------------------------------%
poly_coeff=polyfit(logh,logerr,1);
ordine_convergenza = poly_coeff(1);
fprintf('The'' global order of convergence L2 calcolato con polyfit è: %.3f\n', ordine_convergenza);
log_E_fit = polyval(poly_coeff, logh);

% Per plottare, riportiamo la retta teorica nel mondo "normale" con l'esponenziale
E_fit = exp(log_E_fit);

poly_coeffH=polyfit(logh,logerrH,1);
ordine_convergenzaH = poly_coeffH(1);
fprintf('The'' global order of convergence H2 calcolato con polyfit è: %.3f\n', ordine_convergenzaH);
log_E_fitH = polyval(poly_coeffH, logh);

% Per plottare, riportiamo la retta teorica nel mondo "normale" con l'esponenziale
E_fitH = exp(log_E_fitH);

figure;
% Plot dell'errore effettivo
loglog(h, errors, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;

% Plot della retta trovata con polyfit/polyval
loglog(h, E_fit, 'r--', 'LineWidth', 1.5);

loglog(h, errorsH, 'yo-', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;

% Plot della retta trovata con polyfit/polyval
loglog(h, E_fitH, 'g--', 'LineWidth', 1.5);
legend('Errore L^2', 'Fit O(h^2)', 'Errore H^1', 'Fit O(h^1)', 'Location', 'best');

grid on;
title('Test di Convergenza Errore L^2 & H^1');
xlabel('Dimensione elemento (h)');
ylabel('Errore L^2');

figure;
ur=ureal(coord(:,1),coord(:,2));
trisurf(triang,coord(:,1),coord(:,2),ur-uh);
title('grafico errore')

figure;
trisurf(triang,coord(:,1),coord(:,2),ur);
title('real surface')

figure;
trisurf(triang,coord(:,1),coord(:,2),uh);
title('aprroximated surface')
%
%plot_tri(coord,triang,1,'node_index_l')
%plot_tri(coord,triang,1,'edges')
%plot_tri(coord,triang,1,'face_index')