function [Vflatten, V] = eh()

    V = zeros(50, 50);
    i = linspace(-0.5,0.5,50);
    j = linspace(-0.5,0.5,50);
    
    for m = 1:50
        for n = 1:50
            x = i(m);
            y = j(n);          
            V(m, n) = 0.5*((1-2*y)*((1-2*x)/abs(1-2*y)*asinh(abs((1-2*y)/(1-2*x)))+(1-2*x)/abs(1-2*x)*asinh(abs((1-2*x)/(1-2*y))))+(1+2*y)*((1-2*x)/abs(1+2*y)*asinh(abs((1+2*y)/(1-2*x)))+(1-2*x)/abs(1-2*x)*asinh(abs((1-2*x)/(1+2*y))))+(1-2*y)*((1+2*x)/abs(1-2*y)*asinh(abs((1-2*y)/(1+2*x)))+(1+2*x)/abs(1+2*x)*asinh(abs((1+2*x)/(1-2*y))))+(1+2*y)*((1+2*x)/abs(1+2*y)*asinh(abs((1+2*y)/(1+2*x)))+(1+2*x)/abs(1+2*x)*asinh(abs((1+2*x)/(1+2*y)))));
        end
    end

    surf(meshgrid(i), meshgrid(j).', V);
    xlabel('X');
    ylabel('Y');
    zlabel('V(X, Y)');
    Vflatten = reshape(V, 1, 2500);
    %V = V';
    
end

function val = mrdivide(a, b)
    if b == 0
        val = 0;
    else
        val = builtin('mrdivide', a, b);
    end
end