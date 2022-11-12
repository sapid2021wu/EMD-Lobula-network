% Computational code underlying the main findings in the manuscript 
% "Bioinspired figure-ground discrimination via visual motion smoothing".
% Submitted to Plos Computational Biology.
% Nov 12, 2022.


clear all
close all

% stimulus pattern.  
step = 0.01;                 % in [s].
phi = 0.33;                  % 1 [pixel] = 0.33 [deg]
nx0 = round(90/phi);         % transforming [deg] to [pixel].  
ny0 = round(180/phi); 
bar_width = round(25/phi);   % bar width, transforming [deg] to [pixel]. 
dot_size = 8*1;             % dot size, in [pixel].
%dot_size = 8*3;

% parameters defining types of stimuli, for instance:
% two types of 'Bar on Ground' stimuli (abs(v_b)>=100; abs(v_o)>=100):
% with stationary background: motion_type = 1; v_b = 0; v_o = 200; 
% with moving background: motion_type = 1; v_b = -200; v_o = 200; 
% two types of 'Theta Figure in Moving Background' stimuli (abs(v_b)>=100; abs(v_o)>=100):
% with small dot texture: motion_type = -1; v_b = -200; v_o =-200; dot_size = 8*1;
% with large dot texture: motion_type = -1; v_b = -200; v_o =-200; dot_size = 8*3;

 motion_type = -1;  
 v_b = -200;                  % background velocity in [pixel/s].
 v_o = -200;                  % bar velocity in [pixel/s]. 
 

% reading the stimulus pattern
[patt, bar_location] = stimulus(nx0, ny0, bar_width, dot_size, step, v_b, v_o, motion_type);
input_GT = bar_location;     % groundtruth
% EMD parameters
EMD_nx=size(patt, 1)-1;  
EMD_ny=size(patt, 2)-1;  
frame=size(patt, 3);
% low-pass filtering. 
tauL = 0.050;                 % the filter's time constant in [s].  
dL = tauL/step;                       
a_low(1, 1) = 1/(dL+1);
a_low(1, 2)= 1-a_low(1, 1);
% high-pass filtering. 
tauH = 0.250;                 % the filter's time constant in [s]. 
dH = tauH/step;                     
b_high(1, 1) = dH/(dH+1);
b_high(1, 2) = b_high(1, 1);

% EMD array
Fh=zeros(EMD_nx+1, EMD_ny+1);           
On=zeros(EMD_nx+1, EMD_ny+1); 
Off=zeros(EMD_nx+1, EMD_ny+1);
Fd_On=zeros(EMD_nx+1, EMD_ny+1);      
Fd_Off=zeros(EMD_nx+1, EMD_ny+1);       
H_On = zeros(EMD_nx, EMD_ny);        % output of EMD array in the ON pathway
H_Off = zeros(EMD_nx, EMD_ny);       % output of EMD array in the OFF pathway
threshold = 0.002;                   % output threshold of EMD output. 
 
receptive_s = 5;           
sigma_s = receptive_s/6;           
filter_EMD_s = 1.0*fspecial('gaussian', receptive_s, sigma_s);  % receptive field of Bl and Br units.

% lobula network parameters 
Lc_tau = 5;          % membrane time constant of lobula units in [ms].
beta = 0.5;          % steepness of the activation function of lobula units.
Vhalf0 = -40;        % half-activation voltage of the activation function of lobula units in [mV]. 
Lc_step = 0.4;       % integration time step of the lobula units in [ms].
Con_EMD_L1 = 150;    % synaptic weight from the EMD array to the Br and Bl modules. 
Lc_nx = EMD_nx;      % width of individual lobula modules
Lc_ny = EMD_ny;      % height of individual lobula modules

v = -50+ones(3*Lc_nx, Lc_ny);                     % membrane potentials of the Br, Bl, and Bm modules.
v_2 = -50+zeros(3*Lc_nx, Lc_ny, frame);           % reacording membrane potentials of Br, Bl, and Bm.  
Con_Lc_plot = zeros(3*Lc_nx, Lc_ny, frame);       % output of Br, Bl, and Bm.  
RF_Bm = 3;                               
filter_Bm = fspecial('gaussian', RF_Bm, RF_Bm/6); % receptive field of Bm units.

D_v = -50*ones(3*Lc_nx, Lc_ny);                   % membrane potentials of the Lr, Ll, and Lm modules.
D_v_2 = -50*ones(3*Lc_nx, Lc_ny, frame);          % reacording membrane potentials of Lr, Ll, and Lm.
filter_D_Lc = (-0.05)*[1 0 -1; 1 0 -1; 1 0 -1];   % receptive field of the units in Lr, Ll, and Lm.

Conduct_inputBr = zeros(Lc_nx, Lc_ny, frame);     % input of Br.
EMD_o = zeros(EMD_nx, EMD_ny, frame);             % output of the EMD array

for t = 2:frame
        Fh = b_high(1, 1)*(patt(:, :, t)-patt(:, :, t-1))+b_high(1, 2)*Fh; 
        On = Rect(1, (Fh+0.1*patt(:, :, t)), 0);
        Off = Rect(-1, (Fh+0.1*patt(:, :, t)), 0.05); 
        Fd_On = a_low(1, 1)*On+a_low(1, 2)*Fd_On; 
        Fd_Off = a_low(1, 1)*Off+a_low(1, 2)*Fd_Off;
        % ON pathway 
        H_On(1:(EMD_nx-1), 1:(EMD_ny-1)) = Fd_On(1:(EMD_nx-1), 1:(EMD_ny-1)).*On(1:(EMD_nx-1), 2:EMD_ny)-On(1:(EMD_nx-1), 1:(EMD_ny-1)).*Fd_On(1:(EMD_nx-1), 2:EMD_ny);
        % OFF pathway 
        H_Off(1:(EMD_nx-1), 1:(EMD_ny-1)) = Fd_Off(1:(EMD_nx-1), 1:(EMD_ny-1)).*Off(1:(EMD_nx-1), 2:EMD_ny)-Off(1:(EMD_nx-1), 1:(EMD_ny-1)).*Fd_Off(1:(EMD_nx-1), 2:EMD_ny);
       
        H_On_1 = H_On;   
        H_Off_1 = H_Off;
        
        H_On_1(abs(H_On_1)<=threshold)=0;
        H_Off_1(abs(H_Off_1)<=threshold)=0;
        EMD_o(:, :, t) = H_On_1 + H_Off_1;
                      
        Hon_L = H_On_1;
        Hon_R = Hon_L;
        Hoff_L = H_Off_1;  
        Hoff_R = Hoff_L;
        
        Hon_L(Hon_L>=0)=0;
        Hon_R(Hon_R<0)=0;
        Hoff_L(Hoff_L>=0)=0;
        Hoff_R(Hoff_R<0)=0;
             
        Con_R = conv2(Con_EMD_L1*(Hon_R + Hoff_R), filter_EMD_s, 'same');   % positive components in the input of Br.
        Con_L = conv2(-Con_EMD_L1*(Hon_L + Hoff_L), filter_EMD_s, 'same');  % negative components in the input of Br.
        Conduct_inputBr(:, :, t) = Con_R-Con_L; 
        
        [Con_Lc, v_new, D_v_new] = LobulaUnits(Lc_tau, Lc_step, Lc_nx, v, D_v, Con_R, Con_L, filter_Bm, filter_D_Lc, Vhalf0, beta);
                                     
        v = v_new; 
        D_v = D_v_new; 
                                   
        v_2(:, :, t) = v_new;                
        D_v_2(:, :, t) = D_v_new;           
        Con_Lc_plot(:, :, t) = Con_Lc;           
end

% calculating the F-measures
F_5 = zeros(5, frame-1);
t_m = round(frame/2);
for i = 1:(frame-1)                  
        
     imGT = uint8(255*input_GT(:, :, i));                  % groundtruth
                        
     imBinary_EMD_01 = mat2gray(EMD_o(:, :, i));     
     imBinary_EMD_0 = im2bw(imBinary_EMD_01, 0.5);         % segmentation threshold was set as a 50% maximum.
     imBinary_EMD = uint8(255*imBinary_EMD_0);             % output of the EMD array.
        
     Conduct_R01 = mat2gray(Conduct_inputBr(:,:, i));
     imBinary_Conduct_R0 = im2bw(Conduct_R01, 0.5);
     imBinary_Conduct_R = uint8(255*imBinary_Conduct_R0);  % input of the Br module.
     
     Conduct_L01 = mat2gray(-Conduct_inputBr(:,:, i));             
     imBinary_Conduct_L0 = im2bw(Conduct_L01, 0.5);
     imBinary_Conduct_L = uint8(255*imBinary_Conduct_L0);  % input of the Bl module.
     
     imBinary_L1_00 = mat2gray(Con_Lc_plot(:, :, i+1));      
     imBinary_L1_0 = im2bw(imBinary_L1_00, 0.5); 
     imBinary_L1 = uint8(255*imBinary_L1_0);               % output of the Br, Bl, and Bm modules.
     
     if i==t_m
         imBinary_EMD_m = imBinary_EMD;
         imBinary_Conduct_R_m = imBinary_Conduct_R;
         imBinary_Conduct_L_m = imBinary_Conduct_L;
         imBinary_L1_m = imBinary_L1;         
     end
      
     F_EMD = F_measure(imBinary_EMD, imGT);                       % F-measure at the output of the EMD array
     F_Conduct_R = F_measure(imBinary_Conduct_R, imGT);           % F-measure at the input of the Br module
     F_Conduct_L = F_measure(imBinary_Conduct_L, imGT);           % F-measure at the input of the Bl module
     F_Br = F_measure(imBinary_L1(1:Lc_nx, :), imGT);             % F-measure at the output of the Br module
     F_Bl = F_measure(imBinary_L1((Lc_nx+1):2*Lc_nx, :), imGT);   % F-measure at the output of the Bl module
     
     F_5(:, i) = [F_EMD; F_Conduct_R; F_Conduct_L; F_Br; F_Bl]; 
end
                      
%===============================================
%=======================Activities of the EMD array, Br, Bl, Lr, and Ll========================
figure(2)
set(gcf,'Units','normalized','Position', [0.02 0.1 0.44 0.5]);
%----------------------
subplot(3, 2, 1)
M = moviein(frame);
for j = 1:frame
     imagesc(EMD_o(:, :, j), [-0.23, 0.23]);  % only one frame!
     M(:, j) = getframe;
end
title('Output of the ON+OFF EMD array');
bar1 = colorbar('YTick', -0.23:0.23:0.23,  'YTickLabel', {'-0.23', '', '+0.23'},'FontSize', 8); 
set(bar1, 'position', [0.5020    0.7075    0.02    0.16]);
%----------------------
subplot(3, 2, 3)
M = moviein(frame);
for j = 1:frame
     imagesc(v_2((Lc_nx+1):2*Lc_nx, :, j), [-74, -10]); 
     M(:, j) = getframe;
end
title('Membrane potentials of the Bl module');
%----------------------
subplot(3, 2, 4)
M = moviein(frame);
for j = 1:frame
     imagesc(v_2(1:Lc_nx, :, j), [-74, -10]); 
     M(:, j) = getframe;
end
title('Membrane potentials of the Br module');
bar4 = colorbar('YTick', -74:32:-10,  'YTickLabel', {'-74', '-42', '-10'},'FontSize', 8); 
set(bar4, 'position', [0.9420    0.4071    0.02    0.16]);
%----------------------
subplot(3, 2, 5)
M = moviein(frame);
for j = 1:frame
     imagesc(D_v_2((Lc_nx+1):2*Lc_nx, :, j), [-70, -14] ); 
     M(:, j) = getframe;
end
title('Membrane potentials of the Ll module');
%----------------------
subplot(3, 2, 6)
M = moviein(frame);
for j = 1:frame
     imagesc(D_v_2(1:Lc_nx, :, j), [-70, -14] ); 
     M(:, j) = getframe;
end
title('Membrane potentials of the Lr module');
bar6 = colorbar('YTick', -70:28:-14,  'YTickLabel', {'-70', '-42', '-14'},'FontSize', 8); 
set(bar6, 'position', [0.9420    0.1068    0.02    0.16]);
%===================================================
%=========== Membrane potentials of 4 units in Br ==============
figure(3)
colormap('gray');
set(gcf, 'Units', 'normalized', 'Position', [0.45 0.06 0.3 0.4]);
%----------------------
subplot(2, 1, 1)
imagesc(v_2(1:Lc_nx, :, t_m), [-74, -10]); 
hold on;
plot(15, 23, 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', [1 1 1]);
hold on;
plot(35, 23, 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'color', [0 0 1], 'MarkerFaceColor', [0 0 1], 'MarkerEdgeColor', [1 1 1]);
hold on;
plot(55, 23, 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'color', [0 0.8 0], 'MarkerFaceColor', [0 0.8 0], 'MarkerEdgeColor', [1 1 1]);
hold on;
plot(75, 23, 'o', 'LineWidth', 1.5, 'MarkerSize', 5, 'color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', [1 1 1]);
hold on;
text(25, -6, 'Membrane potentials of the Br module');
text(20, -2, '(when the bar passed the middle of the visual field)');
bar1 = colorbar('YTick', -74:32:-10,  'YTickLabel', {'-74', '-42', '-10'},'FontSize', 8); 
set(bar1, 'position', [0.9420    0.5803    0.02    0.16]);
%----------------------
subplot(2, 1, 2)
Vm_Br(:, 1) = v_2(23, 15, :);
Vm_Br(:, 2) = v_2(23, 35, :);
Vm_Br(:, 3) = v_2(23, 55, :);
Vm_Br(:, 4) = v_2(23, 75, :);

plot(1:frame, Vm_Br(:, 1), '-', 'color', [0.6 0.6 0.6], 'Linewidth', 2);
hold on;
plot(1:frame, Vm_Br(:, 2), '-', 'color', [0 0 1], 'Linewidth', 2);
hold on;
plot(1:frame, Vm_Br(:, 3), '-', 'color', [0 0.8 0], 'Linewidth', 2);
hold on;
plot(1:frame, Vm_Br(:, 4), '-', 'color', [1 0 0], 'Linewidth', 2);
hold on;
text(22, 8, 'Time courses of the membrane potentials of the 4 units marked above');
xlabel('Time (s)', 'Fontsize', 8);  
ylabel('Vm (mV)', 'Fontsize', 8);  
set(gca, 'XLim', [0, 240], 'XTick', 0:80:240, 'XTickLabel', {'0.0', '0.8', '1.6', '2.4'}, 'Fontsize', 8);
set(gca, 'YLim', [-80, 0], 'YTick', -80:40:0, 'YTickLabel', {'-80', '-40', '0'}, 'Fontsize', 8);
%===================================================
%=========== Snapshots of the segmented foreground and background ==============
figure(4)
colormap('gray');
%----------------------
subplot(3, 2, 1)
imagesc(imBinary_EMD_m);    % the middle frame.
text(-30, -10, 'Snapshots of the segmented foreground and background', 'FontSize', 10);
text(-15, -4, '(when the bar passed the middle of the visual field)', 'FontSize', 8);
text(95, 20, 'output of', 'FontSize', 8);
text(95, 27, 'the EMD array', 'FontSize', 8);
set(gca, 'position', [0.3600    0.7093    0.3347    0.2157]);
%----------------------
subplot(3, 2, 3)
imagesc(imBinary_Conduct_L_m); 
hold on;
[fy, fx] = find(input_GT(:, :, t_m), 1, 'first');
[Ly, Lx] = find(input_GT(:, :, t_m), 1, 'last');
plot([fx  fx], [fy Ly], '-r');
hold on;
plot([Lx  Lx], [fy Ly], '-r');
hold on;
plot([fx  Lx], [fy fy], '-r');
hold on;
plot([fx  Lx], [Ly Ly], '-r');
hold on;  
text(30, -5, 'input of Bl', 'FontSize', 8);
axis off;
%----------------------
subplot(3, 2, 4)
imagesc(imBinary_Conduct_R_m); 
hold on;
[fy, fx] = find(input_GT(:, :, t_m), 1, 'first');
[Ly, Lx] = find(input_GT(:, :, t_m), 1, 'last');
plot([fx  fx], [fy Ly], '-r');
hold on;
plot([Lx  Lx], [fy Ly], '-r');
hold on;
plot([fx  Lx], [fy fy], '-r');
hold on;
plot([fx  Lx], [Ly Ly], '-r');
hold on;  
text(30, -5, 'input of Br', 'FontSize', 8);
axis off;
%----------------------
subplot(3, 2, 5)
imagesc(imBinary_L1_m((Lc_nx+1):2*Lc_nx, :)); 
hold on;
[fy, fx] = find(input_GT(:, :, t_m), 1, 'first');
[Ly, Lx] = find(input_GT(:, :, t_m), 1, 'last');
plot([fx  fx], [fy Ly], '-r');
hold on;
plot([Lx  Lx], [fy Ly], '-r');
hold on;
plot([fx  Lx], [fy fy], '-r');
hold on;
plot([fx  Lx], [Ly Ly], '-r');
hold on;  
text(30, -5, 'output of Bl', 'FontSize', 8);
axis off;
%----------------------
subplot(3, 2, 6)
imagesc(imBinary_L1_m(1:Lc_nx, :)); 
hold on;
[fy, fx] = find(input_GT(:, :, t_m), 1, 'first');
[Ly, Lx] = find(input_GT(:, :, t_m), 1, 'last');
plot([fx  fx], [fy Ly], '-r');
hold on;
plot([Lx  Lx], [fy Ly], '-r');
hold on;
plot([fx  Lx], [fy fy], '-r');
hold on;
plot([fx  Lx], [Ly Ly], '-r');
hold on;  
axis off;
text(30, -5, 'output of Br', 'FontSize', 8);
%===================================================
%===================F-measure===================
figure(5)
set(gcf,'Units','normalized','Position', [0.5 0.25 0.4 0.34]);
t = size(F_5, 2);
%----------------------
subplot(1, 2, 1)
F_EMD = F_5(1, :);
F_input_Bl = F_5(3, :);
F_output_Bl = F_5(5, :);

plot(1:t, F_EMD, '-', 'color', [0 0.4 0]);
hold on;
plot(1:t, F_input_Bl, '-', 'color', [0 0.7 1]);
hold on;
plot(1:t, F_output_Bl, '-', 'color', [0 0 0.7]);
hold on;
leg = legend('output of EMD array', 'input of Bl', 'output of Bl');
set(leg, 'position', [0.28, 0.55, 0.1, 0.05]);
text(120, 1.05, 'F-measures throughout the entire stimulus presentation period');
text(15, 0.95, 'the EMD-Bl part');
xlabel('Time (s)','FontSize',10); 
ylabel('F-measure', 'FontSize',10); 
set(gca, 'XLim', [0, 250], 'XTickLabel', {'0.0', '0.5', '1.0', '1.5', '2.0', '2.5'},  'FontSize',8);
set(gca, 'YLim', [0, 1.0]);
%----------------------
subplot(1, 2, 2)
F_EMD = F_5(1, :);
F_input_Br = F_5(2, :);
F_output_Br = F_5(4, :);

plot(1:t, F_EMD, '-', 'color', [0 0.4 0]);
hold on;
plot(1:t, F_input_Br, '-', 'color', [0 0.7 1]);
hold on;
plot(1:t, F_output_Br, '-', 'color', [0 0 0.7]);
hold on;
leg = legend('output of EMD array', 'input of Br', 'output of Br');
set(leg, 'position', [0.65, 0.2, 0.1, 0.05]);
text(15, 0.95, 'the EMD-Br part');
xlabel('Time (s)', 'FontSize',10); 
ylabel('F-measure', 'FontSize',10); 
set(gca, 'XLim', [0, 250], 'XTickLabel', {'0.0', '0.5', '1.0', '1.5', '2.0', '2.5'},  'FontSize',8);
set(gca, 'YLim', [0, 1.0]);
%===================================================
%===================================================
