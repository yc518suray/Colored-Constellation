% domain_coloring.m
% This program demonstrates the idea of domain coloring in constellation points.
%
% Author: Raymond Su, raymondsu0110@gmail.com

%% main

clear
rng();

% ===== constellation generation ===== %
mod_size = 64; % constellation size
Xdim = 32;
Ydim = 32;
Xdim_data = 16;
Ydim_data = 16;
%s = randi(mod_size, [Ydim, Xdim]) - 1;
s = (7 + 7j) * ones(Ydim_data, Xdim_data);
s(1, :) = -(7 + 7j) * ones(1, Xdim_data);
s(Ydim_data, :) = -(7 - 7j) * ones(1, Xdim_data);
%s = qammod(s, mod_size, "gray");
s = [zeros(round((Ydim - Ydim_data) / 2), Xdim); zeros(Ydim_data, round((Xdim - Xdim_data) / 2)), ...
      s, zeros(Ydim_data, round((Xdim - Xdim_data) / 2)); zeros(round((Ydim - Ydim_data) / 2), Xdim)];

% ===== presentation of Tx symbols ===== %
show_color_constellation(s, mod_size, "Tx DD-domain matrix", 1);
x = ifft(s, Xdim, 2) * sqrt(Xdim);

% ===== channel ===== %
xt = x(:);
h = [1, 0.28 - 0.15j, -0.35 + 0.32j, 0, 0, 0.15 + 0.18j];
nu = [0, 6, -4, 0, 0, 2];
rt = zeros(1, length(xt) + length(h));
for k = 1: length(h)
    rt = rt + h(k) * exp(1j * 2 * pi * nu(k) * (1: Xdim * Ydim + length(h)) / (Xdim * Ydim)) .* ...
          [zeros(1, k - 1), xt.', zeros(1, length(h) - k + 1)];
end
rt = rt(1: Xdim * Ydim);
r = reshape(rt, [Ydim, Xdim]);
y = fft(r, Xdim, 2) / sqrt(Xdim);

% ===== presentation of Rx symbols ===== %
show_color_constellation(y, mod_size, "Rx DD-domain matrix", 0);

clear


%% function definitions
function show_color_constellation(s, N, name, opt)
    % plot the colored constellation
    % s -> matrix of symbols
    % N -> constellation size
    % name -> name of the matrix to be plotted
    % opt -> option to enable plot of constellation diagram, 1 for YES and 0 for NO
    %
    % Note: the constellation should be square
    
    % constellation settings
    if opt == 1
        constella = gen_constellation(N);
        constella_ext = [constella, zeros(length(constella(:, 1)), 1)];
        constella_ext = [constella_ext; zeros(1, length(constella_ext(1, :)))];
        [Ny, Nx] = size(constella_ext);
        [xc, yc] = meshgrid(0: Nx - 1, 0: Ny - 1);
        
        colored_constella = mat2rgb(constella_ext);
    end
    
    % symbol matrix settings
    s = [s, zeros(length(s(:, 1)), 1)];
    s = [s; zeros(1, length(s(1, :)))];
    [My, Mx] = size(s);
    [x, y] = meshgrid(0: Mx - 1, 0: My - 1);
    
    colored_s = mat2rgb(s);
    
    % plot the constellation
    if opt == 1
        figure();
        c = pcolor(xc, yc, zeros(size(constella_ext)));
        set(c, "CData", colored_constella);
        axis equal; axis([0, Nx - 1, 0, Ny - 1]);
        xticks([0.5: Nx - 0.5]); yticks([0.5: Ny - 0.5]);
        xticklabels(string(-Nx + 2: 2: Nx - 2));
        yticklabels(string(-Ny + 2: 2: Ny - 2));
        xlabel("real part"); ylabel("imaginary part");
    
        title_str = sprintf("Constellation, size = %d", N);
        title(title_str);
    end
    
    % plot the matrix of symbols
    figure();
    c = pcolor(x, y, zeros(size(s)));
    set(c, "CData", colored_s, "EdgeColor", "None");
    axis equal; axis([0, Mx - 1, 0, My - 1]);
    xticks([0.5: nextpow2(Mx) - 1: Mx - 0.5]); yticks([0.5: nextpow2(My) - 1: My - 0.5]);
    xticklabels(string(0: nextpow2(Mx) - 1: Mx - 1));
    yticklabels(string(0: nextpow2(My) - 1: My - 1));
    xx = x(1: My - 1, 1: Mx - 1);
    yy = y(1: My - 1, 1: Mx - 1);
    xx = reshape(xx, [1, (Mx - 1) * (My - 1)]) + 0.3;
    yy = reshape(yy, [1, (Mx - 1) * (My - 1)]) + 0.5;
    
    % optional (constellation point annotations)
    %text(xx, yy, string(s(1: end - 1, 1: end - 1)), "FontSize", 24 / sqrt(0.6 * (My + Mx - 2)));
    
    % title
    title_str = sprintf("%s, constellation size = %d", name, N);
    title(title_str);
end

function rgb = mat2rgb(X)
    % convert matrix X to rgb matrix

    persistent MAX_SATURATION;
    if isempty(MAX_SATURATION)
        MAX_SATURATION = 0;
    end
    
    h = (angle(X) + pi) / (2 * pi); % hue (color)
    s = abs(X);
    if MAX_SATURATION == 0
        MAX_SATURATION = max(s(:));
    end
    s = s / MAX_SATURATION; % saturation (grayness)
    v = ones(size(X)); % value (brightness)
    
    rgb = hsv2rgb(cat(3, h, s, v));
end

function C = gen_constellation(N)
    % generate a typical constellation diagram
    % N -> constellation size
    
    K = sqrt(N);
    x = -K + 1: 2: K - 1;
    C = ones(K, 1) * x + 1j * x' * ones(1, K);
end