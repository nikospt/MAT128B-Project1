% MAT 128B: Project 1
% UC Davis Winter 2020
% Nikos Trembois, Caitlin Brown, and Shuai Zhi

clc; close all; clearvars

global c xRange yRange pts bsave
% Constants that will be used for z^2 - c plots
c = [0.36 + 0.1i, -.123 - .745i,-.749, -1.25];
xRange = [1, 1, 1, 1.6]; % Range of x values for plotting window
yRange = [1, 1, 1, 0.7]; % Range of y values for plotting window
pts = 100; % Number of points in x and y directions of plot
bsave = 0; % boolean value (0,1) to save plots to file or not

%% Part 1: Fractals
phi = @(z,c) z^2; % Defining the function
% Iterating the Julia Set with function defined at bottom of script
M = FilledJuliaSet(phi,1,1,pts,0,'nocolor',2);

figure(); hold on
%title('Filled Julia Set of $\phi = z^2$','Fontsize',16,'Interpreter','Latex')
xlabel('\Re','Fontsize',18); ylabel('\Im','Fontsize',18)
colormap([1 0 0; 1 1 1]); 
image( [-1 1], [-1 1], M')
axis xy; axis equal; axis([ -1 1 -1 1])
hold off
if bsave == 1
    saveas(gcf,'../Figures/UnitDisk.png')
end

%% Part 2:  Fractals
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2 + c;
M = cell(length(c),1);
for i = 1:length(c)
    M{i} = FilledJuliaSet(phi,xRange(1),yRange(1),pts,c(i),'nocolor',2);
end

% Plotting the maps of orbits
for i = 1:length(c)
    figure(); hold on
    % Generate the title string based on constant value
    if (real(c(i)) > 0 )
        stitle = strcat('Filled Julia Set of $z^2 + ',num2str(c(i)),'$');
    else
        stitle = strcat('Filled Julia Set of $z^2 ',num2str(c(i)),'$');
    end
    %title(stitle,'Interpreter','Latex','FontSize',24) % Create Title
    colormap([1 0 0; 1 1 1]); % Define colors for map
    image( [-xRange(i) xRange(i)], [-yRange(i) yRange(i)], M{i}') % create map
    axis xy; axis equal; ax = gca;
    ax.XLim = [-xRange(i) xRange(i)]; ax.YLim = [-yRange(i) yRange(i)];
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18); ylabel('\Im','Fontsize',18)
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    hold off
    if bsave == 1 % Save plot to file if bsave = 1
        ssave = strcat('../Figures/FilledJulia',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end

%% Part 3: Julia Sets
clearvars -except c xRange yRange pts bsave
psi = @(z,c) sqrt(z - c); % Define plot 
x = zeros(pts,length(c)); % initialize a x vector for each constant
y = zeros(pts,length(c)); % initialize a y vector for each constant
for k = 1:length(c)
    clear z;
    z = c(k);
    for j = 1:1000 % Use 1,000 iterations
        x(j,k) = real(z(j));
        y(j,k) = imag(z(j));
        % Randomly choose positive or negative root with equal weighting
        if randi([0 1],1,1) == 1
            z(j+1) = psi(z(j),c(k));
        else
            z(j+1) = -psi(z(j),c(k));
        end
    end
end

% Plot the Julia sets for different constant values
for i = 1:length(c)
    figure(); hold on
    % Generate the title string based on constant value
    if (real(c(i)) > 0 )
        stitle = strcat('Julia Set of $z^2 + ',num2str(c(i)),'$');
    else
        stitle = strcat('Julia Set of $z^2 ',num2str(c(i)),'$');
    end
    %title(stitle,'Interpreter','Latex','FontSize',24) % Create title
    scatter(x(:,i),y(:,i),'filled') % Plot the Julia Set
    % Plot Formatting
    ax = gca;
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18); ylabel('\Im','Fontsize',18)
    hold off
    if bsave == 1 % Save plot to file if bsave = 1
        ssave = strcat('../Figures/Julia',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end

%% Fractal Dimension
% Using square domain to capture fractal dimension more easily
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2 + c;
range = 2; % define a range used for both x and y, for square spacing
M = cell(length(c),1);
for i = 1:length(c)
    M{i} = FilledJuliaSet(phi,range,range,200,c(i),'nocolor',2);
end
% Calculating the fractal dimension
for i = 1:length(c)
    r = 2*range/pts; % calculate result
    fprintf('For c = (%4.2f, %4.2fi)',real(c(i)),imag(c(i)))
    FractalDimension(M{i},r)
end

%% Part 5: Connectivity of the Julia Set
clearvars -except c xRange yRange pts bsave
psi = @(z) z^2 + 3; % Calculate inverse iteration undder some function psi
z = 0; % Initialize value
max_iter = 1000; % Use 1000 iterations to determine if orbit diverges
for i = 1:max_iter
    z = psi(z); % Calculate the orbit
    if abs(z) > 100 % Considered divergent when magnitude is greater than 100
        fprintf('The orbit diverged after %i iterations, the set is not connected\n',i)
        break
    end
end
if abs(z) < 100
    fprintf('The set did not diverge after %i iterations\n',max_iter)
    fprintf('It is reasonable to assume the Julia set is connected\n')
end
    
%% Part 6: Coloring Divergent Orbits
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2+ c;
for i = 1:length(c)
    % Once again using the FilledJuliaSet function is used with options 
    % 'colored' to color plots based on iterations to diverge
    % past an absolute value of 100
    M{i} = FilledJuliaSet(phi, xRange(i), yRange(i), 200, c(i), 'colored', 100);
end

% Plotting the Diverging orbits
for i = 1:length(c)
    figure(); hold on
    % plot map of diverging orbits
    image( [-xRange(i) xRange(i)], [-yRange(i) yRange(i)], M{i}')
    % Plot formatting below
    axis xy; axis equal; ax = gca;
    ax.XLim = [-xRange(i) xRange(i)]; ax.YLim = [-yRange(i) yRange(i)];
    plot(ax.XLim,[0,0],'LineStyle','--','Color',[.5,.5,.5])
    plot([0,0],ax.YLim,'LineStyle','--','Color',[.5,.5,.5])
    xlabel('\Re','Fontsize',18); ylabel('\Im','Fontsize',18)
    colormap(jet(max(max(M{i})))); colorbar; hold off
    if bsave == 1 % Save plot to file if bsave = 1
        ssave = strcat('../Figures/ColoredJulia',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end

%% Part 7: Newton's Method in the Complex Plane
clearvars -except c xRange yRange pts bsave
% Defining function that represents f/f' for a general f = z^n - 1
g = @(z,n) (z^n - 1)/(n*z^(n-1));
a = linspace(-5,5,500); b = linspace(-5,5,500); % set plottig window
nmax = 5; % The maximum polynomial order to plot to
M = cell(nmax-1,1);
for n = 2:nmax % start with order 2 as order 1 is not interesting
    M{n-1} = 100*ones(length(a),length(b)); % initialize values to
    gn = @(z) (z^n - 1)/(n*z^(n-1)); % Redefine g for specific value of n
    for r = 1:length(a)
        for i = 1:length(b)
            z = a(r) + 1i*b(i);
            for j = 1:100
                if abs(z^n-1) > 0.001 % Testing for convergence
                    z = z - gn(z); % Newton's Iteration
                else
                    if j < 3 % Printing values of possible roots
                        fprintf('For n = %i, The root is near %2.4f, %2.4f\n',n,a(r),b(i));
                    end
                    M{n-1}(r,i) = j; % Set value to number of iterations
                    break;
                end
            end
        end
    end
end

% Plotting the Newton iterations
stitle1 = 'Fixed points of';
for i = 1:nmax-1
    figure(); hold on
    % Creating title string
    stitle2 = strcat( ' $z^',num2str(i+1),' - 1$');
    % title(strcat(stitle1,stitle2),'Interpreter','Latex')
    % Plotting map of iterations of convergence
    image( [min(a) max(a)], [min(b) max(b)], M{i}')
    % Plot formatting
    colormap(jet(max(max(M{i})))); colorbar
    axis xy; axis equal; ax = gca;
    ax.XLim = [min(a) max(a)]; ax.YLim = [min(b) max(b)];
    xlabel('\Re','Fontsize',18); ylabel('\Im','Fontsize',18)
    hold off
    if bsave == 1 % Save plot to file if bsave = 1
        ssave = strcat('../Figures/Newton',num2str(i),'.png');
        saveas(gcf,ssave)
    end
end
    

%% Part 8: Mandelbrot Set
clearvars -except c xRange yRange pts bsave
phi = @(z,c) z^2 + c; % function to plot
a = linspace(-1,1,500); % Plotting Window
b = linspace(-1,1,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        z = 0; % start with z = 0
        cm = a(r) + b(i)*1i;
        for j = 1:100
            z = phi(z,cm);
            if abs(z) > 100
                M(r,i) = j;
                break;
            end
        end
    end
end

% Plot Mandelbrot set
figure(); hold on
% title('Mandelbrot Set')
xlabel('Real'); ylabel('Imaginary')
colormap(jet(100)); colorbar
image( [-1 1], [-1 1], M')
axis xy; axis equal; axis([ -1 1 -1 1])
if bsave == 1 % Save plot to file if bsave = 1
    saveas(gcf,'../Figures/Madelbrot.png')
end
%% Part 8 (cont) Zoom in
% zoom in on a fractal by changing limits
clear M;
phi = @(z,c) z^2 + c;
a = linspace(-.3,0.05,500);
b = linspace(0.6,1,500);
M = ones(length(a),length(b));

for r = 1:length(a)
    for i = 1:length(b)
        z = 0;
        cm = a(r) + b(i)*1i;
        for j = 1:100
            z = phi(z,cm);
            if abs(z) > 100
                M(r,i) = j;
                break;
            end
        end
    end
end

figure(); hold on
%title('Mandelbrot Set')
xlabel('Real'); ylabel('Imaginary')
colormap(jet(100)); colorbar
image( [min(a) max(a)], [min(b) max(b)], M')
axis xy; axis('equal')

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
    fprintf(' the fractal dimension is %5.4f\n', d);
end

% This function determines if a point is in a Julia Set
function [ M ] = FilledJuliaSet(phi, xrange, yrange, pts, c, colored, maxvalue)
    a = linspace(-xrange, xrange, pts);
    b = linspace(-yrange, yrange, pts);
    if (strcmp(colored,'colored'))
        M = zeros(length(a), length(b));
    else
        M = ones(length(a), length(b));
    end
    for r = 1:length(a)
        for i = 1:length(b)
            clear z;
            z = a(r) + 1i*b(i);
            for j = 1:100
                z = phi( z, c );
                if abs( z ) > maxvalue
                    if (strcmp(colored,'colored'))
                        % If the colored option is 'colored', then
                        % color the julia set by the number of 
                        % iterations to converge
                        M(r,i) = j;
                    else % otherwise give discrete value
                        M(r,i) = 2;
                    end
                    break
                end
            end
        end
    end
    if (strcmp(colored,'colored'))
        M(M==0) = max(max(M)) + 1;
    end
end