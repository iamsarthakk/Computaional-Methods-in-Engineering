function NotAKnotCubicSpline(x, y, t )
x = sort(x);

[n,~] = size(x);
order = 4*(n-1);
order_linear = 2*(n-1);
A = zeros(order);
B = zeros(order,1);

B(1) = y(1);
for i = 2:2:(order_linear-2)
    B(i) = y((i/2)+1);
    B(i+1) = B(i);
end
B(order_linear) = y(end);
A(1,1) = 1;
A(1,2) = x(1);
A(1,3) = A(1,2)^2;
A(1,4) = A(1,2)^3;
for i = 2:2:(order_linear-2)
    num = 2*i;
    A(i,num-3) = 1;
    A(i,num-2) = x((i/2)+1);
    A(i,num-1) = A(i,num-2)^2;
    A(i,num) = A(i,num-2)^3;
    A(i+1,num+1:num+4) = A(i,num-3:num);
end
i = i+2;
num = 2*i;
A(i,num-3) = 1;
A(i,num-2) = x(end);
A(i,num-1) = A(i,num-2)^2;
A(i,num) = A(i,num-1)^3;
i = i+1;
for i = 1:(n-2)
    num = i + order_linear;
    A(num,(4*i)-2) = 1;
    A(num,(4*i)-1) = 2*x(i+1);
    A(num,(4*i)) = 3*((x(i+1))^2);
    A(num,((4*i)+2):((4*i)+4)) = -1*A(num,((4*i)-2):(4*i));
end
for i = 1:(n-2)
    num = num+1;
    A(num,(4*i)-1) = 2;
    A(num,(4*i)) = 6*x(i+1);
    A(num,((4*i)+3):((4*i)+4)) = -1*A(num,((4*i)-1):(4*i));
end
A(end-1,4) = 6;
A(end-1,8) = -6;
A(end,end-4) = 6;
A(end,end) = -6;
CP = A\B;

[N,~] = size(x);
Y = zeros(N,1);
CC = zeros(n-1,2);
CC(1:end,1) = CP(1:4:end);
CC(1:end,2) = CP(2:4:end);
CC(1:end,3) = CP(3:4:end);
CC(1:end,4) = CP(4:4:end);


for i=1:length(t)
    for j=1:length(x)
        if t(i)<=x(j)
            break
        end
    end
    k=j-1;
    a=CC(k,4);
    b=CC(k,3);
    c=CC(k,2);
    d=CC(k,1);
    y_out(i)= a*(t(i))^3 + b*(t(i))^2 + c*(t(i)) + d;
end

% Plotting results:
hold on;
for i = 1:n-1
    curve = @(x) CC(i,1) + CC(i,2).*(x) + CC(i,3).*(x).^2 + CC(i,4).*(x).^3;
    fplot(curve,[x(i),x(i+1)],'y');
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
title('Not-a-Knot Cubic Spline')

%Output
fileID=fopen('output.txt','w');
fprintf(fileID,'%s \n','Not-a-Knot Cubic spline');
for i=1:length(t)
    fprintf(fileID,'%6.4f %6.4f \n',t(i),y_out(i));
end