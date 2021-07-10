%Bessel function is a mathematical function needed to generate the Kaiser window in FIR filter
function [ I_x ] = getBessel(x)% k is considered as 500
summation=1;
for k=1:500
    summation=summation+((1/factorial(k))*(x/2).^k).^2;
end
I_x=summation;
end
