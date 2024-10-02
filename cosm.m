function outputArg = cosm(inputArg)
%COSM cos value of a matrix
%   cosm(A) = (exp(iA) + exp(-iA)) / 2
outputArg = (expm(1i * inputArg) + expm(-1i * inputArg)) / 2;
end