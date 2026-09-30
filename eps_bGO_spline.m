function eps = eps_bGO_spline(w,axis) %d representa el 'nodo' en que quieres saber la permitividad
    % data from Gonzalo
    switch axis
            case 'X'                      % 100              
                eps = eps_A_bGO_3parOsc(w);  
            case 'Z'                      % 010             
                eps = eps_B_bGO_3parOsc(w);  
            case 'Y'                      % 001             
                eps = eps_C_bGO_3parOsc(w); 
    end
end
function 
data = readmatrix('bGO_undoped_Jose.txt');
eps_xx = interp1(data(:,1),data(:,2),w,'spline') + 1i*interp1(data(:,1),data(:,3),w,'spline');
eps_xy = interp1(data(:,1),data(:,4),w,'spline') + 1i*interp1(data(:,1),data(:,5),w,'spline');
eps_yy = interp1(data(:,1),data(:,10),w,'spline') + 1i*interp1(data(:,1),data(:,11),w,'spline');
eps_zz = interp1(data(:,1),data(:,18),w,'spline') + 1i*interp1(data(:,1),data(:,19),w,'spline');
% epsT_bGO = zeros(3,3,Nw);
epsT(1,:,:) = [eps_xx;eps_xy;zeros_w(w);]; 
epsT(2,:,:) = [eps_xy;eps_yy;zeros_w(w);]; 
epsT(3,:,:) = [zeros_w(w);zeros_w(w);eps_zz;];


end