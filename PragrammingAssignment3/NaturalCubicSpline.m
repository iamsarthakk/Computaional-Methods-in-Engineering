function NaturalCubicSpline(x,y,t)

n = length(x);
 

h = zeros(n-1,1);
for j = 1:n-1
 h(j) = x(j+1) - x(j);
end
 

A = zeros(n);
 

A(1,1)= 1;
A(n,n) = 1;
 
for i = 2:n-1
 A(i,i-1) = h(i-1);
 A(i,i) = 2*(h(i-1)+h(i));
 A(i,i+1) = h(i);
end
 

b = zeros(n,1);
 
for i = 2:n-1
 b(i) = (3/h(i))*(y(i+1)-y(i)) - (3/h(i-1))*(y(i)-y(i-1));
end
 

cj = A\b;
 
bj = zeros(n-1,1);
for i = 1:n-1
 bj(i) = (1/h(i))*(y(i+1)-y(i)) - (1/3*h(i))*(2*cj(i)+cj(i+1));
end
 

dj = zeros(n-1,1);
for i = 1:n-1
 dj(i) = (1/(3*h(i))) * (cj(i+1)-cj(i));
end
 

% Given a value on the x-axis, inter, that we wish to know the y-value of,
% we must first find in which interval inter is. We will use bisection
% search to accomplish this. Interval must be ascending or descending.
 
for i=1:length(t)
    for j=1:length(x)
        if t(i)<=x(j)
            break
        end
    end
    k=j-1;
    a=dj(k);
    b=cj(k);
    c=bj(k);
    d=y(k);
    y_out(i)= a*(t(i)-x(k))^3 + b*(t(i)-x(k))^2 + c*(t(i)-x(k)) + d;
end
 
% Plotting results:
hold on;
for i = 1:n-1
    curve = @(u) y(i) + bj(i).*(u-x(i)) + cj(i).*(u-x(i)).^2 + dj(i).*(u-x(i)).^3;
    fplot(curve,[x(i),x(i+1)],'m');
    hndl=get(gca,'Children');
    set(hndl,'LineWidth',2);
    hold on
end
scatter(x,y,50,'r','filled')
grid on;
xlim([min(x) max(x)]);
ylim([min(y)-1 max(y)+1]);
xlabel('x');
ylabel('y');
title('Natural Cubic Spline')

%Output
fileID=fopen('output.txt','w');
fprintf(fileID,'%s \n','Natural Cubic spline');
for i=1:length(t)
    fprintf(fileID,'%6.4f %6.4f \n',t(i),y_out(i));
end