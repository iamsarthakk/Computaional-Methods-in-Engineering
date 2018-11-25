roots = [10 20 30 40 50];
n = 5;
close_guess = 32;
t1 = 20;
    for i_pr = 1:n
        t2 = abs(roots(i_pr) - close_guess);
        if(t2<t1)
            t1 = t2;
            r_req = roots(i_pr);
            disp(t2);
        end
    end
    disp(r_req);