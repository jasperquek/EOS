% ---- Initializing online parameters -------
    G= 10^11;%shear modulus Stromboli = 10^9.4; Gelatin = 5000
    mu= 10^3.4; %viscosity Stromboli 10^3
    rho= 2400; %magma density Stromboli = 2600
    rc =20; %conduit radius Stromboli = 3m
    mass0 = 10.^(6.7);%5,5.5,5.7,6; %kg Stromboli = 780 for analogue experiment
    
        online_data = [G, mu, rho, rc, mass0];
        
save('Online_variables','online_data')