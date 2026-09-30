function eps = eps_CrSBr_spline(w,Axis) %d representa el 'nodo' en que quieres saber la permitividad
	
    switch Axis
            case 'X'                      % 100, x           
                eps = eps_A_CrSBr_spline(w);  
            case 'Y'                      % 010, y          
                eps = eps_B_CrSBr_spline(w);  
            case 'Z'                      % 001, z            
                eps = eps_C_CrSBr_spline(w); 
    end
end

function eps = eps_A_CrSBr_spline(w)
dataReal = readmatrix('eps1_RT_aaxis.dat');
dataImag = readmatrix('eps2_RT_aaxis.dat');

eps=1;  %auxiliar para salida de función además del plot

epsr = interp1(dataReal(:,1),dataReal(:,2),w,'spline');
epsi = interp1(dataImag(:,1),dataImag(:,2),w,'spline');
eps = epsr+1i*epsi;
% disp(eps);
end

function eps = eps_B_CrSBr_spline(w)
dataReal = readmatrix('eps1_RT_baxis.dat');
dataImag = readmatrix('eps2_RT_baxis.dat');

eps=1;  %auxiliar para salida de función además del plot

epsr = interp1(dataReal(:,1),dataReal(:,2),w,'spline');
epsi = interp1(dataImag(:,1),dataImag(:,2),w,'spline');
eps = epsr+1i*epsi;
% disp(eps);
end

function eps = eps_C_CrSBr_spline(w)
dataReal = readmatrix('eps1_RT_caxis.dat');
dataImag = readmatrix('eps2_RT_caxis.dat');

eps=1;  %auxiliar para salida de función además del plot

epsr = interp1(dataReal(:,1),dataReal(:,2),w,'spline');
epsi = interp1(dataImag(:,1),dataImag(:,2),w,'spline');
eps = epsr+1i*epsi;
% disp(eps);
end