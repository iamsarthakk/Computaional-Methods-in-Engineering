
fileid = fopen('c:\\t.txt', 'r');
I = fscanf(fileid, '%f', inf);
fclose(fileid);
n = I(1);
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
        tmp = sqrt((tmp(1))^2 - (y1(iter))^2);
        iter = iter + 1;
    end
    y_tr = trans_pose(y,n,1);
    %disp(y);
    lambda1 = (mul(mul(y_tr , A,1,n,n), y,1,n,1))/(mul(y_tr, y,1,n,1));
    e = abs((lambda1 - lambda)/(lambda));
    fprintf('Lamda\t%.6f %d %.9f\n',lambda1,i,e);
    lambda = lambda1;
    y = z / abs(tmp);
    %disp(lambda);
    %disp(i);
    i = i+1;
end
end
%disp(y);

end
    
    
