function Problem1()
x = 1;

y25_5iter = series(x,.25,5);
y50_5iter = series(x,.50,5);
y75_5iter = series(x,.75,5);

y25_3iter = series(x,.25,3);
y50_3iter = series(x,.50,3);
y75_3iter = series(x,.75,3);

disp("First 5")
disp([y25_5iter,y50_5iter,y75_5iter])

disp("***********")

disp("First 3")
disp([y25_3iter,y50_3iter,y75_3iter])

end


function temp = series(x,y,iter)
L = 2;
w = 1;
t1 = 50;
t2 = 150;

count = 1;
n = 1;
theta_out = 0;
while count <= iter
    theta = ( (((-1)^(n+1)+1) / n) * sin(n*pi*x/L) * ( sinh(n*pi*y/L) / sinh(n*pi*w/L) ) );
    if theta ~= 0
        count = count + 1;
        theta_out = theta_out + theta;
    end
    n = n + 1;
end
theta_out = (2/pi) * theta_out;

temp = theta_out*(t2-t1) + t1;
end


