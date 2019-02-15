function Problem1()
x = 1;

y25_5iter = series(x,.25,500000000000000);
y50_5iter = series(x,.50,5);
y75_5iter = series(x,.75,5);

y25_3iter = series(x,.25,3);
y50_3iter = series(x,.50,3);
y75_3iter = series(x,.75,3);

disp([y25_5iter,y50_5iter,y75_5iter])
disp([y25_3iter,y50_3iter,y75_3iter])

end


function temp = series(x,y,iter)
L = 2;
w = 1;
t1 = 50;
t2 = 150;

theta_ = 0;
for i = 1:iter
    theta = (2/pi) * ((-1^(iter+1)+1)*sin(iter*pi*x/L)*(sinh(iter*pi*y/L)/sinh(iter*pi*w/L)));
    if theta ~= 0
        disp('Non zero')
    end
    theta_ = theta_ + theta;
end

temp = theta_*(t2-t1) + t1;
end