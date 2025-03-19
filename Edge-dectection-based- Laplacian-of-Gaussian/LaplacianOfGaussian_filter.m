%  Laplacian of Gaussian filter
clc;clear;
datain = load('Synthetic data.txt');
nx=200; ny=200;
data = datain(:,3);
x0 = datain(:,1);
y0 = datain(:,2);
xx = reshape(x0,nx,ny);
yy = reshape(y0,nx,ny);
data0 = reshape(data,nx,ny);

k=1;
kernelSize = 3;
sigma=0.01;
[m,n] = size(data0);
e1 = zeros(m,n);
rr = 2:m-1; cc=2:n-1;

% Creating Nomrmalize the laplacian of gaussian kernel
LoGKernel = fspecial( 'log', kernelSize, sigma );

% Filter the data to generate a response to the laplacian of gaussian
filteredData = imfilter(data0, LoGKernel,'same', 'replicate','conv');

b = filteredData;
thresh = 2400;
% Look for the zero crossings:  +-, -+ and their transposes

e2 = zeros(m,n);
[rx,cx] = find( b(rr,cc) < 0 & b(rr,cc+1) > 0 ...
    & abs( b(rr,cc)-b(rr,cc+1) ) > thresh );   % [- +]
e2((rx+1) + cx*m) = 1;

[rx,cx] = find( b(rr,cc-1) > 0 & b(rr,cc) < 0 ...
    & abs( b(rr,cc-1)-b(rr,cc) ) > thresh );   % [+ -]
e2((rx+1) + cx*m) = 1;

[rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc) > 0 ...
    & abs( b(rr,cc)-b(rr+1,cc) ) > thresh);   % [- +]'
e2((rx+1) + cx*m) = 1;

[rx,cx] = find( b(rr-1,cc) > 0 & b(rr,cc) < 0 ...
    & abs( b(rr-1,cc)-b(rr,cc) ) > thresh);   % [+ -]'
e2((rx+1) + cx*m) = 1;

[rx,cx] = find( b(rr-1,cc-1) > 0 & b(rr,cc) < 0 ...
    & abs( b(rr-1,cc-1)-b(rr,cc) ) > thresh);   % [+ -]'
e2((rx+1) + cx*m) = 1;

[rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc+1) > 0 ...
    & abs( b(rr,cc)-b(rr+1,cc+1) ) > thresh);   % [+ -]'
e2((rx+1) + cx*m) = 1;

% Create the Laplacian of Gaussian kernel with different kernelSize
for i=1:9
LoGKernel = fspecial( 'log', kernelSize, sigma );
filteredDate = imfilter(data0, LoGKernel,'same', 'replicate','conv');

[m,n] = size(data0);
e = zeros(m,n);
rr = 2:m-1; cc=2:n-1;
b = filteredDate;
[rx,cx] = find( b(rr,cc) < 0 & b(rr,cc+1) > 0 ...
    & abs( b(rr,cc)-b(rr,cc+1) ) > thresh );   % [- +]
e((rx+1) + cx*m) = 1;
[rx,cx] = find( b(rr,cc-1) > 0 & b(rr,cc) < 0 ...
    & abs( b(rr,cc-1)-b(rr,cc) ) > thresh );   % [+ -]
e((rx+1) + cx*m) = 1;
[rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc) > 0 ...
    & abs( b(rr,cc)-b(rr+1,cc) ) > thresh);   % [- +]'
e((rx+1) + cx*m) = 1;
[rx,cx] = find( b(rr-1,cc) > 0 & b(rr,cc) < 0 ...
    & abs( b(rr-1,cc)-b(rr,cc) ) > thresh);   % [+ -]'
e((rx+1) + cx*m) = 1;
[rx,cx] = find( b(rr,cc) < 0 & b(rr+1,cc+1) > 0 ...
    & abs( b(rr,cc)-b(rr+1,cc+1) ) > thresh);   % [+ -]'
e((rx+1) + cx*m) = 1;
[rx,cx] = find( b(rr-1,cc-1) > 0 & b(rr,cc) < 0 ...
    & abs( b(rr-1,cc-1)-b(rr,cc) ) > thresh);   % [+ -]'
e((rx+1) + cx*m) = 1;

outdata=e(:);
% selecting the coordinates x, y 
sel_x=[];sel_y=[];sel_data=[];
new_sol_i = 0; 
ii=abs(m*n);
 for i=1:m*n
        if  outdata(i)>0
        new_sol_i=new_sol_i + 1;
        sel_x(new_sol_i,1) = x0(i,1); sel_y(new_sol_i,1) = y0(i,1); 
        sel_data(new_sol_i,1) = outdata(i,1);
    end;
  end;

figure(1);
subplot(3,3,k)
scatter(sel_x,sel_y,0.1,'MarkerEdgeColor','k','MarkerFaceColor','k')
axis equal
axis square
title(sprintf('Kernel size: %2.0f',kernelSize));
xlabel({'X(km)'})
ylabel({'Y(km)'})
set(gca,'FontName','Times New Roman', 'FontSize',12)
set(get(gca,'XLabel'),'FontSize',12,'FontName','Times New Roman')
set(get(gca,'YLabel'),'FontSize',12,'FontName','Times New Roman')
xlim([0 200])
ylim([0 200])
box on
k=k+1;
kernelSize =kernelSize+2;
 end
exportgraphics(gcf,'model2_diffkernel.jpg','Resolution',300)


