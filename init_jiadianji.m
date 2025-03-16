function pop = init_jiadianji(pop_size,dim,ub,lb)
p = zeros(pop_size,dim);
prime_number_min = dim*2+3;
while 1
    if isprime(prime_number_min) == 1
        break;
    else
        prime_number_min=prime_number_min+1;
    end
end
for i =1:pop_size
    for j = 1:dim
        r = mod(2*cos(2*pi*j/prime_number_min)*i,1);
        p(i,j) = lb+r*(ub-lb);
    end
end
pop = p;
end