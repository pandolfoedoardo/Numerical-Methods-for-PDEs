function err=evalErr(Nelem,coord,triang,Area,uh,ureal)
err=0;
for i = 1:Nelem
    elemerr=0; %local error
        for j = 1:3 %trapeizodal rule
           nj=triang(i,j);
           coordj=coord(nj,:);
           elemerr=elemerr + (Area(i)/3)*(uh(nj)-ureal(coordj(1),coordj(2)))^2; 
        end
err=err+elemerr;  %global error
end
err=sqrt(err);
end