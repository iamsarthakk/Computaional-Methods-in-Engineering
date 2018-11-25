function[] = qrd()
[filename,filepath] = uigetfile('.txt','File Selector');
s= fullfile(filepath,filename);
fileid = fopen(s , 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I(1,1);
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
        Q = qr_dec(A,n);
        %disp(Q);
        R = qr_dec1(A,Q,n);
        disp(R);
end
function[q] = qr_dec(a,n)
q = zeros(n,n,'double');
q(:,1) = a(:,1)/no_rm(a(:,1),n);
for k = 1:(n-1)
    z1 = double.empty(n,0);
    sum = 0;
    for i = 1:k
        sum = sum + (mul(trans_pose(a(:,k+1),n,1),q(:,i),1,n,1))*q(:,i);
    end
    z1 = a(:,k+1) - sum;
    q(:,k+1) = z1/no_rm(z1,n);
end
end
function[r] = qr_dec1(a,q,n)
r = zeros(n,n,'double');
q_t = trans_pose(q,n,n);
r = mul(q_t,a,n,n,n);
end
function[tmp] = no_rm(y1,n)
 tmp = y1(1);
    iter = 2;
    while(iter<=n)
        tmp = sqrt((tmp(1))^2 + (y1(iter))^2);
        iter = iter + 1;
    end
end
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