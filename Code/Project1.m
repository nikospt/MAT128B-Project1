% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and Shuai Zhi

%% Part 1: Fractals
phi = @(z) z^2;
a = linspace(-1,1,500);
b = linspace(-1,1,500);
M = ones(length(a),length(b));
for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = a(r) + 1i*b(i);
        for j = 1:100
            z(j+1) = phi(z(j));
            if abs(z(j+1)) > 2
                M(r,i) = 2;
                break;
            end
        end
    end
end

figure(); hold on
title('Filled Julia Set of $\phi = z^2$','Fontsize',16,'Interpreter','Latex')
xlabel('\Re','Fontsize',18)
ylabel('\Im','Fontsize',18)
colormap([1 0 0; 1 1 1]);
image( [-1 1], [-1 1], M)
axis xy
axis('equal')
axis([ -1 1 -1 1])
hold off
saveas(gcf,'../Figures/UnitDisk.png')

%% Part 2:  Fractals
clear M
phi = @(z,c) z^2 + c;
rl = -1.6; ru = -rl;
il = -.7;  iu = -il;
a = linspace(rl,ru,500);
b = linspace(il,iu,500);
%c = [0.36 + 0.1i, -.123 - .745i,-.749,-.25+.25i, -1.25]; % 
c = -1.25;
M = cell(length(c),1);
for k = 1:length(c)
    M{k} = ones(length(a),length(b));
    for r = 1:length(a)
        for i = 1:length(b)
            clear z;
            z = a(r) + 1i*b(i);
            for j = 1:100
                z(j+1) = phi(z(j),c(k));
                if abs(z(j+1)) > 2
                    M{k}(r,i) = 2;
                    break;
                end
            end
        end
    end
end


FractalDimension(M{end},.02)

for i = 1:length(c)
    figure(); hold on
    colormap([1 0 0; 1 1 1]);
    image( [rl ru], [il iu], M{i})
    axis xy
    hold off
end

%% Part 3: Julia Sets
psi = @(z,c) sqrt(z - c);
x = zeros(100,4);
y = zeros(100,4);
c = [0.36 + 0.1i, -.123 - .745i,-.749,-.25+.25i, -1.25];
for k = 1:length(c)
    clear z;
    z = c(k);
    for j = 1:10000
        x(j,k) = real(z(j));
        y(j,k) = imag(z(j));
        if randi([0 1],1,1) == 1
            z(j+1) = psi(z(j),c(k));
        else
            z(j+1) = -psi(z(j),c(k));
        end
    end
end

for i = 1:length(c)
    figure(); hold on
    title('Julia Set of $z^2 + c$','Interpreter','Latex','FontSize',24)
    xlabel('Real')
    ylabel('Imaginary')
%     ax = gca;
%     ax.XAxisLocation = 'origin';
%     ax.YAxisLocation = 'origin';
    scatter(x(:,i),y(:,i),'filled')
    hold off
end

%% Part 4: Computing the Fractal Dimension

%% Part 5: Connectivity of the Julia Set

%% Part 6: Coloring Divergent Orbits

%% Part 7: Newton's Method in the Complex Plane

%% Part 8: Mandelbrot Set

phi = @(z,c) z^2 + c;
a = linspace(-1,1,500);
b = linspace(-1,1,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = 0;
        c = a(r) + b(i)*1i;
        for j = 1:100
            z(j+1) = phi(z(j),c);
            if abs(z(j+1)) > 100
                M(r,i) = j;
                break;
            end
        end
    end
end

figure(); hold on
title('Mandelbrot Set')
xlabel('Real')
ylabel('Imaginary')
colormap(jet(100))
image( [-1 1], [-1 1], M)
colorbar
axis xy
axis('equal')
axis([ -1 1 -1 1])
%% Part 8 (cont) Zoom in
% zoom in on a fractal by changing limits
% Weird it seems to change significantly when using higher resolution in
% this section from what is shown in the previous plot
phi = @(z,c) z^2 + c;
a = linspace(-1,-.8,500);
b = linspace(0,-.2,500);
clear M;
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        clear z;
        z = 0;
        c = a(r) + b(i)*1i;
        for j = 1:100
            z(j+1) = phi(z(j),c);
            if abs(z(j+1)) > 2
                M(r,i) = j;
                break;
            end
        end
    end
end

figure(); hold on
title('Mandelbrot Set')
xlabel('Real')
ylabel('Imaginary')
colormap(jet(100))
image( [min(a) max(a)], [min(b) max(b)], M)
colorbar
axis xy
axis('equal')

%% Functions
function [result] = R(x,y)
    result = sqrt(x^2 + y^2);
end

function [result] = T(x,y)
    if x > 0
        t = atan(y/x);
    elseif x == 0
        if y > 0
            t = pi/2;
        else
            t = 3*pi/2;
        end
    else
        t = atan(y/x) + pi;
    end
    result = t;
end

function FractalDimension(M,r)
    n = length(M);
    N = 0; % Initialize value of N for box-counting
    for i = 1:n
        for j = 1:n
            if M(i,j) == 1
                N = N+1;
            end
        end
    end
    d = log(N)/log(1/r);
    fprintf('The fractal dimension is %5.4f\n', d);
end
