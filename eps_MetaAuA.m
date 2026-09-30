function c = eps_MetaAuA(w) %d representa el 'nodo' en que quieres saber la permitividad
M = readmatrix('PermittivitiesMetaAu.dat');

c = 1;  %auxiliar para salida de función además del plot
freq = 1./(M(:,1)*10^(-4));
epsreal = M(:,2);
epsimag = M(:,3);

epsr = interp1(freq,epsreal,w,'spline');
epsi = interp1(freq,epsimag,w,'spline');
c = epsr+1i*epsi;
end