function sol = eval_anal_Dani( x, z, g )

% ---------------------------------------------------------------------------
% ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION:
%
% BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
% FAR FIELD FLOW - VISCOSITIES - GEOMETRY
%
% Yolanda Deubelbeiss,       10.03.2007
% ---------------------------------------------------------------------------

% INPUT:
gr  = 0;                        % Simple shear: gr=1, er=0
er  = -1;                       % Strain rate
mm  =  1;                       % Viscosity of matrix
mc  =  1e3;                     % Viscosity of clast
rc  =  0.5*g.xmax;              % Radius of clast

A   =   mm.*(mc-mm)./(mc+mm);
i   =   sqrt(-1);


% --------------------------------------------------------------
% PRESSURE CALCULATION OUTSIDE OF AN INCLUSION IN THE Z-PLANE
% --------------------------------------------------------------

% INSIDE CLAST
if sqrt(x^2 + z^2)<=rc
    
    Z       =   x + i*z;
    sol.P   =   0;   % if you want you can add NaN, but according to Schmid's thesis it's zero inside
    
    % VELOCITY
    V_tot          =  (mm/(mc+mm))*(i*gr+2*er)*conj(Z)-(i/2)*gr*Z;
    sol.vx         =  real(V_tot);
    sol.vz         =  imag(V_tot);
    
    
    % OUTSIDE CLAST, RESP. MATRIX
else
    Z              =   x + i*z;
    % PRESSURE
    sol.P          =   -2.*mm.*(mc-mm)./(mc+mm).*real(rc^2./Z.^2.*(i*gr+2*er));
    
    % VELOCITY
    phi_z          = -(i/2)*mm*gr*Z-(i*gr+2*er)*A*rc^2*Z^(-1);
    d_phi_z        = -(i/2)*mm*gr + (i*gr+2*er)*A*rc^2/Z^2;
    conj_d_phi_z   = conj(d_phi_z);
    psi_z          = (i*gr-2*er)*mm*Z-(i*gr+2*er)*A*rc^4*Z^(-3);
    conj_psi_z     = conj(psi_z);
    
    V_tot          = (phi_z- Z*conj_d_phi_z - conj_psi_z) / (2*mm);
    sol.vx         =  real(V_tot);
    sol.vz         =  imag(V_tot);
    
end


end

