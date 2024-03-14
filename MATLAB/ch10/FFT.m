function g = FFT(f)
    N = length(f);
    if N == 1
        g = f;
        return;
    end

    k = log2(N);
    %不满2的幂次，补齐
    if ceil(k)~=k
        n = 2^ceil(k);
        tmp = zeros(n-N,1);
        f = [f;tmp];
    else
        n = N;
    end

    f1 = f(1:2:end);
    f2 = f(2:2:end);
    g1 = FFT(f1);
    g2 = FFT(f2);
    
    w_n = exp(-1i*2*pi/n);
    w = 1;
    g = zeros(size(f));
    for k = 1:n/2
        g(k) = (g1(k)+w*g2(k))/2;
        g(k+n/2) = (g1(k)-w*g2(k))/2;
        w = w*w_n;
    end

    if n~=N
        g = g(1:N);
    end
end