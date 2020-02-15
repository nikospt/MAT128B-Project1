% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and Shuai Zhi

%% Part 1: Fractals
clearvars
phi = @(z,c) z^2;
M = FilledJuliaSet(phi,1,1,100,0);

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
clearvars
phi = @(z,c) z^2 + c;
xrange = 2; yrange = 2;
c = [0.36 + 0.1i, -.123 - .745i,-.749,-.25+.25i, -1.25];
M = cell(length(c),1);
for i = 1:length(c)
    M{i} = FilledJuliaSet(phi,xrange,yrange,100,c(i));
end

for i = 1:length(c)
    figure(); hold on
    colormap([1 0 0; 1 1 1]);
    image( [-xrange xrange], [-yrange yrange], M{i})
    axis xy
    axis equal
    ax = gca;
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    hold off
end

%% Part 3: Julia Sets
clearvars
psi = @(z,c) sqrt(z - c);
x = zeros(100,4);
y = zeros(100,4);
c = [0.36 + 0.1i, -.123 - .745i,-.749,-.25+.25i, -1.25];
for k = 1:length(c)
    clear z;
    z = c(k);
    for j = 1:1000
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
    scatter(x(:,i),y(:,i),'filled')
    ax = gca;
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    hold off
end

%% Part 4: Computing the Fractal Dimension

% FractalDimension(M{end},.02)

%% Part 5: Connectivity of the Julia Set
clearvars
psi = @(z) z^2 + 3;
z = 0;
max_iter = 1000;
for i = 1:max_iter
    z = psi(z);
    if abs(z) > 100
        fprintf('The orbit diverged after %i iterations, the set is not connected\n',i)
        break
    end
end
if abs(z) < 100
    fprintf('The set did not diverge after %i iterations\n',max_iter)
    fprintf('It is reasonable to assume the Julia set is connected\n')
end
    

%% Part 6: Coloring Divergent Orbits
clearvars
phi = @(z,c) z^2 - c;
rl = -1.6; ru = -rl;
il = -.7;  iu = -il;
a = linspace(rl,ru,100);
b = linspace(il,iu,100);
c = [0.36 + 0.1i, -.123 - .745i,-.749,-.25+.25i, -1.25]; 
%c = -1.25;
M = cell(length(c),1);
for k = 1:length(c)
    M{k} = ones(length(a),length(b));
    for r = 1:length(a)
        for i = 1:length(b)
            clear z;
            z = a(r) + 1i*b(i);
            for j = 1:100
                z(j+1) = phi(z(j),c(k));
                if abs(z(j+1)) > 100
                    M{k}(r,i) = j;
                    break;
                end
            end
        end
    end
end

for i = 1:length(c)
    figure(); hold on
    image( [rl ru], [il iu], M{i})
    axis xy
    ax = gca;
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    colormap(jet(max(max(M{i}))))
    colorbar
    hold off
end

%% Part 7: Newton's Method in the Complex Plane
clearvars
g = @(z,n) (z^n - 1)/(n*z^(n-1));
a = linspace(-2,2,500);
b = linspace(-2,2,500);
nmax = 3;
M = cell(nmax,1);
for n = 1:3
    M{n} = 100*ones(length(a),length(b));
    gn = @(z) (z^n-1)/(n*z^(n-1));
    for r = 1:length(a)
        for i = 1:length(b)
            z = a(r) + 1i*b(i);
            for j = 1:100
                diff = z^n - 1;
                if abs(z^n-1) > 0.001
                    z = z - gn(z);
                else
                    M{n}(r,i) = j;
                    break;
                end
            end
        end
    end
end

for i = 1:nmax
    figure(); hold on
    image( [min(a) max(a)], [min(b) max(b)], M{i})
    colormap(jet(max(max(M{i}))))
    colorbar
    axis xy
    axis equal
    xlabel('\Re','Fontsize',18)
    ylabel('\Im','Fontsize',18)
    hold off
end
    

%% Part 8: Mandelbrot Set
clearvars
phi = @(z,c) z^2 + c;
a = linspace(-1,1,100);
b = linspace(-1,1,100);
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
a = linspace(-1,-.8,100);
b = linspace(0,-.2,100);
clear M;
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

function [ M ] = FilledJuliaSet(phi, xrange, yrange, pts, c)
    a = linspace(-xrange, xrange, pts);
    b = linspace(-yrange, yrange, pts);
    for k = 1:length(c)
        M = ones(length(a), length(b));
        for r = 1:length(a)
            for i = 1:length(b)
                clear z;
                z = a(r) + 1i*b(i);
                for j = 1:100
                    z(j+1) = phi( z(j), c(k) );
                    if abs( z(j+1) ) > 2
                        M(r,i) = 2;
                        break;
                    end
                end
            end
        end
    end
end