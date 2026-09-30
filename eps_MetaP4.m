function eps = eps_MetaP4(w,Axis) 
eps = 1; 
switch Axis
        case 'X'                                  
            load('Pattern4_effEpsilon.mat','wn','eps_X');
            eps = interp1(wn,eps_X,w,'spline');
        case 'Y'                                  
            load('Pattern4_effEpsilon.mat','wn','eps_Y');
            eps = interp1(wn,eps_Y,w,'spline');
        case 'Z1'                                 
            load('Pattern4_effEpsilon.mat','wn','eps_Z1');
            eps = interp1(wn,eps_Z1,w,'spline');
        case 'Z2'                                 
            load('Pattern4_effEpsilon.mat','wn','eps_Z2');
            eps = interp1(wn,eps_Z2,w,'spline');
end
end