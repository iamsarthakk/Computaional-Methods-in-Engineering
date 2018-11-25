function[] = power_method
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s , 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I(1);
m = size(I);
max_iter = I(m(1) - 2);
rel_error = I(m(1) - 1);
close_guess = I(m);
 A = zeros(n,n,'double');
        i = 1;
        while(i<=n)
            j = 1;
            while(j<=n)
                A(i,j) = I((i-1)*(n-1) + (i+j));
                j = j + 1;
            end
            i = i + 1;
        end      
z = double.empty(n,0);
i1 = 2;
z(1) = 1;
while(i1<= n)
    z(i1) = 0;
    i1 = i1 + 1;
end
e = 100;
y1 = double.empty(n,0);
y = trans_pose(z,1,n);
%disp(z_tr);
lambda = 0;
i = 1;
while((e>rel_error)&&(i<=max_iter))
    z = A * y;
    y1 = trans_pose(z,n,1);
    tmp = y1(1);
    iter = 2;
    while(iter<=n)
        tmp = sqrt((tmp(1))^2 + (y1(iter))^2);
        iter = iter + 1;
    end
    y_tr = trans_pose(y,n,1);
    %disp(y);
    lambda1 = (mul(mul(y_tr , A,1,n,n), y,1,n,1))/(mul(y_tr, y,1,n,1));
    e = abs((lambda1 - lambda)/(lambda1));
    fprintf('Lamda\t%.6f %d %.9f\n',lambda1,i,e);
    lambda = lambda1;
    y = z / abs(tmp);
    %disp(lambda);
    %disp(i);
    i = i+1;
end
end
%disp(y);
function[c] = mul(a,b,m,n,p)
c = zeros(m,p,'double');
for i = 1:m
    for j = 1:p
        sum = 0;
        for k = 1:n
            sum = sum + a(i,k)*b(k,j);
        end
        c(i,j) = sum;
    end
end
end
function[c] = trans_pose(a,m,n)
c = zeros(n,m,'double');
for i = 1:m
    for j = 1:n
            c(j,i) = a(i,j);
        
    end
end
end
    
    
