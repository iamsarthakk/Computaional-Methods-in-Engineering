[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s, 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I(1) - 1;
z = double.empty(n+1,0);
y = double.empty(n+1,0);
m = I(2*n + 4);
z_u = double.empty(m,0);
y_u = double.empty(m,0);
j = 2;
i = 3;
while(j<=(2*n+2))
    z(j/2) = I(j);
    j = j + 2;
end
while(i<=(2*n + 3))
    y((i-1)/2) = I(i);
    i = i + 2;
end
tmp1 = 2*n + 5;
for k = 1:m
   z_u(k) = I(tmp1);
   tmp1 = tmp1 + 1;
end
syms x;
syms sum ;
syms f(x);
sum = 0;
for I = 1:(n+1)  
    L_i = 1;
    for J = 1:(n+1)
        if(J~=I)
             L_i = L_i*(x - z(J))/(double(z(I) - z(J)));
        end
    end
    sum = sum + double(y(I))*L_i;
end
f(x) = sum;
%disp(f(x));
for m5 = 1:m
    y_u(m5) = double(f(z_u(m5)));
end
[filename,filepath] = uigetfile('.txt','File Selector');
s1= fullfile(filepath,filename);
fileid1 = fopen(s1 , 'w');
fprintf(fileid1,'Lagrange \r\n');
fprintf(fileid1,'X\tf(X) \r\n');
formatspec1 = '%0.4f\t%0.4f\r\n';
for ind = 1:m
fprintf(fileid1,formatspec1,z_u(ind),y_u(ind));
end
x = -1.5:0.1:1.5;
y = f(x);
figure
plot(x,y,'g')
hold on
for tmp2 = 1:n+1
    plot(z(tmp2),f(z(tmp2)),'r*')
    hold on
end