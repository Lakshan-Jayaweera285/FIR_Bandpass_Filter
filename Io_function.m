function [outputa]=Io_function(x)
outputa=0
for i=1:10
    outputa=outputa+(1/factorial(i)*(x/2).^i).^2
end
outputa=outputa+1

end

