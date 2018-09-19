// From https://fr.wikipedia.org/wiki/Algorithme_de_Remez
a = 0
b = 1 - %eps

format('v',21)
dispall=%f

for n = 2:10
    
    if dispall
        printf("N = %d\n",n)
    end
    //////////////////////////////////////////////:
    xk = zeros(n+2,1)
    
    for k = 1:(n+2)
        xk(k) = ((a+b)/2) + ((a-b)/2) * cos(((k-1)*%pi)/(n+1))
    end
    
    if dispall
        disp(xk)
    end
    
    //////////////////////////////////////////////:
    M = zeros(n+2, n+2)
    
    for i = 1:(n+1)
        M(:,i) = xk.^(i-1)
    end
    for i = 1:(n+2)
        M(i, n+2) = (-1).^i
    end
    
    if dispall
        disp(M)
    end
    
    //////////////////////////////////////////////:
    func = zeros(n+2, 1)
    for i = 1:(n+2)
        //func(i) = exp(xk(i))
        func(i) = 1 + xk(i) - 2.^xk(i)
    end
    
    if dispall
        disp(func)
    end
    
    //////////////////////////////////////////////:
    [x0,kerA]=linsolve(M,func)
    
    //////////////////////////////////////////////:
    if dispall
        disp(x0)
    end
    
    for i = 1:(n+1)        
        printf("inline constexpr static double GetCoefficient%d_%d() {\n", n+1, i-1);
        printf("    return %.20e;\n", -x0(i));
        printf("}\n");
    end
    printf("\n");

end
