function c = eps_ScN_A(w) %d representa el 'nodo' en que quieres saber la permitividad
M = readmatrix('Perm_ScN_A.txt');

c = 1;  %auxiliar para salida de función además del plot
freq = M(:,1);
epsreal = M(:,2);
epsimag = M(:,3);

epsr = interp1(freq,epsreal,w,'spline');
epsi = interp1(freq,epsimag,w,'spline');
c = epsr+1i*epsi;
end