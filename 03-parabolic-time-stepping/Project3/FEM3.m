%%% Finite element code for 2D piecewise linear Galekrin 
%%% -div(Diff grad u)(x)=f(x) with Dirichlet BC 
clear all; close all; clc;
Meshes = {'mesh0','mesh1','mesh2'};
errors= zeros(length(Meshes),1);
h = zeros(length(Meshes),1);
p=0;
for j = 1:length(Meshes)

    % Define the equation BC's and coefficients for mesh3 only
    meshdir = Meshes{j};
    
    [Nelem,Nodes,triang,coord,NDir,DirNod, ureal,force,u_0] = input_data(meshdir);
    Diff = 1;
    % calculate element area and elemental coefficients of basis functions
    % (b,c)
    [Aloc,Bloc,Cloc,Area] = localBasis(Nelem,triang,coord);
    
    %build stiffness matrix (without BCs)
    stiffMat = stiffBuild(Nelem,Nodes,triang,Bloc,Cloc,Area,Diff);
    MassMat = MassBuild(Nelem, Nodes,triang,Area);
    T=1;
    h(j)=min(2*sqrt(Area));
    dt=0.01*h(j)^2; %first with h
    % Calcolo di quanti passi "interi" servono per coprire T_max
    NMAX = ceil(T / dt); 
    
    % Ricalcolo il dt effettivo per farlo "incastrare" perfettamente
    dt = T/ NMAX;
    
    time= 0; %inizializza il tempo
    uh = zeros(Nodes,NMAX+1); %storia del vettore soluzione
    X = coord(:,1);
    Y= coord(:,2);
    uh(:,1)= u_0(X,Y);
    theta = 0; 
    
    %left hand side
    SysMat = (stiffMat*dt*theta+MassMat);

    errors(j) = 0;

    
 
    
    for i = 1:NMAX
    
        rhs0= (MassMat-(1-theta)*dt*stiffMat)*uh(:,i);
    
        t_old = time;
        f_old = @(x,y) force(x,y,t_old); %fix old time 
        b_old = forcing2d(Nelem, Nodes, triang, coord, Area, f_old);
        
    
        t_new = time + dt;  
        f_new = @(x,y) force(x,y,t_new);
        b_new = forcing2d(Nelem, Nodes, triang, coord, Area, f_new);
    
        rhs1 = dt*(theta*b_new+(1-theta)*b_old);
    
        rhs = rhs0 + rhs1;
        DirVal = zeros(NDir, 1); 
        [SysMat,rhs] = imposeDirPen(NDir,DirNod,DirVal,SysMat,rhs);
    
        uh(:,i+1) = SysMat\rhs;
        time = t_new;
        
        urealt = @(x,y) ureal(x,y,t_new);
        %calcolo dell'errore in norma h^2
        errors(j) = errors(j) + dt*evalErr(Nelem,coord,triang,Area,uh(:,i+1),urealt)^2;
    end
    errors(j) = sqrt(errors(j));
    logh(j)=log(h(j));
    logerr(j)=log(errors(j));
    ur=ureal(coord(:,1),coord(:,2),time);
        p=p+1;
        subplot(3,3,p)
        trisurf(triang,coord(:,1),coord(:,2),ur-uh(:,NMAX+1));
        title('grafico errore')

        p=p+1;
        subplot(3,3,p)
        trisurf(triang,coord(:,1),coord(:,2),ur);
        title('real surface')

        p=p+1;
        subplot(3,3,p)
        trisurf(triang,coord(:,1),coord(:,2),uh(:,NMAX+1));
        title('aprroximated surface')
        fprintf('il numero di iteraziioni per h = %.1f è %.1f \n',h(j),NMAX)
        

end

poly_coeff=polyfit(logh,logerr,1);
ordine_convergenza = poly_coeff(1);
fprintf('The'' global order of convergence L2 calcolato con polyfit è: %.3f\n', ordine_convergenza);
log_E_fit = polyval(poly_coeff, logh);

% Per plottare, riportiamo la retta teorica nel mondo "normale" con l'esponenziale
E_fit = exp(log_E_fit);


figure;
% Plot dell'errore effettivo
loglog(h, errors, 'bo-', 'LineWidth', 2, 'MarkerFaceColor', 'b'); hold on;

% Plot della retta trovata con polyfit/polyval
loglog(h, E_fit, 'r--', 'LineWidth', 1.5);

grid on;
title('Test di Convergenza Errore L^2 - Theta = 0');
xlabel('Dimensione elemento (h)');
ylabel('Errore L^2');


