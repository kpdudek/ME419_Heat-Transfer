function Problem2()

syms theta n w L x y

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% theta = (  ( ((-1)^(n+1)+1) / n ) * sin((n*pi*x)/L) * ( sinh((n*pi*y)/L) / sinh((n*pi*w)/L) )  );
% 
% theta_p = diff(theta,y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% y = 0;
% L = 2;
% w = 1;
% t1 = 50;
% t2 = 150;

% q_prime_y = (pi*cosh((pi*n*y)/L)*sin((pi*n*x)/L)*((-1)^(n + 1) + 1))/(L*sinh((pi*n*w)/L));
% 
% flux = int(q_prime_y,x,0,2)
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Evaluated from x=0 --> x=2
%(2*sin((pi*n)/L)^2*cosh((pi*n*y)/L) - 2*(-1)^n*sin((pi*n)/L)^2*cosh((pi*n*y)/L))/(n*sinh((pi*n*w)/L))

% Evaluated from x=0 --> x=2 and y=0
%q_out_y = (((-1)^n - 1)*(cos(pi*n) - 1))/(n*sinh((pi*n)/2));


flux = series(5);
disp('Flux Over the Lower Surface')
disp(flux)

end



function flux = series(iter)
% L = 2;
% w = 1;
t1 = 50;
t2 = 150;
k = 50;

count = 1;
n = 1;
q_out = 0;
while count <= iter
    q_out_y = (((-1)^n - 1)*(cos(pi*n) - 1))/(n*sinh((pi*n)/2));
    if q_out_y ~= 0
        count = count + 1;
        q_out = q_out + q_out_y;
        %disp(q_out_y)
    end
    n = n + 1;
end
q_out = (2/pi) * q_out;

flux = -k * (q_out*(t2-t1));
end
