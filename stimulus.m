% stimulus pattern

function  [patt, GT] = stimulus(nx0, ny0, bar_width, dot_size, step, v_b, v_o, motion_type)
% step: temporal step of the EMD array.
% dot_size: the size of the randomly distributed black and white dots of the pattern.
% bar_width: bar width. 
% motion_type: 1 (Fourier motion); -1 (theta motion)
% v_b: background velocity in [pixel/sec].
% v_o: bar object velocity in [pixel/sec]. 

% pattern luminance 
% I0: the average luminance.
% I1: the luminance difference.
% Imax = I0+I1     
% Imin = I0-I1     
% the Michelson contrast = I1/I0. 
I0 = 0.5;
I1 = 0.4; 

% size of the stimulus pattern. 
nx = floor(nx0/dot_size);
ny = floor(ny0/dot_size);

rng('shuffle');
Bg0 = randi([0 1], nx, ny);     
Bg0 = Bg0+(Bg0-1);              % elements become -1 or 1
% the foreground bar,  
Fg0 = randi([0, 1], nx, 5*ny);  
Fg0 = Fg0+(Fg0-1);
% the background
Bg0 = I0+I1*Bg0;
Fg0 = I0+I1*Fg0;

Bg = kron(Bg0, ones(dot_size));
Fg = kron(Fg0, ones(dot_size));
v_e = v_o*(motion_type-1);      % element velocity within the theta bar. 

if abs(v_o)<0.00001             % static object.
      delta_t = 2;              % in [s]
else
      delta_t = 1*(size(Bg, 2)-bar_width-5)/abs(v_o);            
end
frame = round(delta_t/step);  
patt_0 = zeros(size(Bg, 1), size(Bg, 2), frame);
ele_m = floor(abs(v_e)*step);     % moving distance of bar elements relative to the bar
patt_0(:, :, 1) = Bg;    

x0 = round(size(Fg, 2)/2);

if v_o<0
      bar_y = size(patt_0, 2)-bar_width;
      patt_0(:, bar_y:(bar_y+bar_width-1), 1) = Fg(:, (x0 - bar_width +1):x0); 
else
      bar_y = 1;     
      patt_0(:, bar_y:(bar_y+bar_width-1), 1) = Fg(:, x0:(x0 + bar_width -1)); % Insert the new foreground
end
bar_m = floor(v_o*step);          % moving distance of the bar 
bg_m = floor(v_b*step);           % moving distance of the background.
Bg1 = zeros(size(Bg));            % the new background

% the location of moving bar
GT_0  = zeros(size(patt_0, 1), size(patt_0, 2), frame);

for t = 2:frame
          % the background,  
          if v_b<0
                  Bg1(:,  1:(size(Bg, 2)-abs(bg_m))) = Bg(:, (abs(bg_m)+1):end);
                  Bg1(:,  (size(Bg, 2)-abs(bg_m)+1):end) = Bg(:, 1:abs(bg_m));
          elseif  abs(v_b)<0.00001        % static background. 
                  Bg1(:, :) = Bg;
          else                        
                  Bg1(:,  1:bg_m) = Bg(:, (size(Bg, 2)-bg_m+1):end);
                  Bg1(:,  (bg_m+1):end) = Bg(:, 1:(size(Bg, 2)-bg_m));
          end   
          patt_0(:, :, t) = Bg1;
          
          % the foreground
           bar_y = bar_y+bar_m;
           
           if v_e<0
                    x0 = x0 + ele_m;
           elseif abs(v_e)<0.00001
                    x0 = x0;
           else
                    x0 = x0 - ele_m;
           end
                      
           patt_0(:, bar_y:(bar_y + bar_width-1), t) = Fg(:, x0:(x0 + bar_width -1));                  
           GT_0(:, bar_y:(bar_y + bar_width-1), t) = 1;
      
           Bg = Bg1;
end

space_time = zeros(frame, size(patt_0, 2));
% obtaining the space-time plot.
for j=1:frame       
         space_time(j, :) = patt_0(1, :, j);
end

Fsize = 12;      
sigma = 3.5;                                
h = fspecial('gaussian', Fsize, sigma);    % 2D Gaussian filter.
K = 6;                                     % in [pixel], for downsampling.

% preprocessing the visual inputs. 
for j=1:frame
          % blurring input images to simulate the fly's optics
          patt_blur = imfilter(patt_0(:, :, j), h);
          patt(:, :, j) =patt_blur(1:K:end, 1:K:end);
end

% getting the bar location.
for j=1:frame
     GT1 = GT_0(1:K:end, 1:K:end, j);      
     GT(:, :, j) = GT1(1:(end-1), 1:(end-1));     % to align with EMD's output.
end

%================stimulus pattern==================
figure(1)
colormap('gray');
set(gcf,'Units','normalized','Position', [0.02 0.7 0.4 0.15]);
%----------------------
subplot(1, 2, 1)
M = moviein(frame);
for j = 1:frame
     imagesc(patt_0(:, :, j), [0, 1]); 
     M(:, j) = getframe;
end
title('Visual stimulus');
%-------------
subplot(1, 2, 2)
imagesc(space_time);
title(['Space-time plot, frame=', num2str(frame)]);
%==================================

      
      
      
      
      
      
      
      
