% numerical integration for the lobula network in one EMD step (i.e., 10 ms)

function  [Con_Lc, v_new, D_v_new] = LobulaUnits(Lc_tau, Lc_step, Lc_nx, v, D_v, Con_R, Con_L, filter_RL, filter_D_Lc, Vhalf0, beta)
% Lc_tau: membrane time constant of lobula units in [ms].             
% step: temporal step of the EMD array in [ms]. 
% Lc_step: integration step of the lobula units in [ms].
% Con_Lc: output of Br, Bl, and Bm. 
% v, v_new: membrane potentials of the Br, Bl, and Bm modules.  
% D_v, D_v_new: membrane potentials of the Lr, Ll, and Lm modules.  
% beta:   steepness of the activation function.
% Vhalf0: half-activation voltage of the activation function in [mV].

E_e = 0;              % reversal potential of excitatory synapses in [mV].
E_i = -80;            % reversal potential of inhibitory synapses in [mV].
Con_L1_L2 = 20;       % synaptic weight from Br, Bl, and Bm to Lr, Ll, and Lm.

i = 1;
for t_Lc=1:25                      
    
    Con_Lc = v;          
    Con_Lc(Con_Lc>=-50) = 1./(1+exp((Vhalf0-Con_Lc(Con_Lc>=-50))/beta)); 
    Con_Lc(Con_Lc<-50) = 0;    
         
    Con_RL = conv2((Con_Lc(1:Lc_nx, :) + Con_Lc((Lc_nx+1):2*Lc_nx, :)), filter_RL, 'same');  
         
    Con_D_Lc(1:Lc_nx, :) = conv2(Con_Lc(1:Lc_nx, :), filter_D_Lc, 'same');
    Con_D_Lc((Lc_nx+1):2*Lc_nx, :) = conv2(Con_Lc((Lc_nx+1):2*Lc_nx, :), filter_D_Lc, 'same');
    Con_D_Lc((2*Lc_nx+1):3*Lc_nx, :) = conv2(Con_Lc((2*Lc_nx+1):3*Lc_nx, :), filter_D_Lc, 'same');
                              
    Con_L1_to_L2_exc = Con_D_Lc;
    Con_L1_to_L2_inh = Con_D_Lc;
    Con_L1_to_L2_exc(Con_L1_to_L2_exc<0) = 0; 
    Con_L1_to_L2_inh(Con_L1_to_L2_inh>0) = 0;                                
                                      
    I_1(1:Lc_nx, :) = Con_R.*(E_e-v(1:Lc_nx, :)) + Con_L.*(E_i-v(1:Lc_nx, :));                                % input to Br.
    I_1((Lc_nx+1):2*Lc_nx, :) = Con_L.*(E_e-v((Lc_nx+1):2*Lc_nx, :)) + Con_R.*(E_i-v((Lc_nx+1):2*Lc_nx, :));  % input to Bl.
    I_1((2*Lc_nx+1):3*Lc_nx, :) = Con_RL.*(E_e-v((2*Lc_nx+1):3*Lc_nx, :));                                    % input to Bm.
    I_2 = Con_L1_L2*(Con_L1_to_L2_exc.*(E_e-D_v) - Con_L1_to_L2_inh.*(E_i-D_v));                              % input to Lr, Ll, and Lm.
    
    v_new = RK4_Lc(Lc_step, v, I_1, Lc_tau);
    D_v_new = RK4_Lc(Lc_step, D_v, I_2, Lc_tau);
  
    v = v_new;
    D_v = D_v_new;
                                                   
end
 
end
