function outputArg = pauli(inputArg)
%PAULI Pauli matrices in quantum mechanics.
%   The input argument can be 0, 1, 2, 3 (numerical values),
%   or '0', 'x', 'y', 'z' (character values). 
%   Returns sigma_0 = [1, 0; 0, 1], 
%   sigma_1 = sigma_x = [0, 1; 1, 0],
%   sigma_2 = sigma_y = [0, -1i; 1i, 0],
%   sigma_3 = sigma_z = [1, 0; 0, -1].
switch inputArg
    case 0
        outputArg = [1, 0; 0, 1];
    case '0'
        outputArg = [1, 0; 0, 1];
    case 1
        outputArg = [0, 1; 1, 0];
    case 'x'
        outputArg = [0, 1; 1, 0];
    case 2
        outputArg = [0, -1i; 1i, 0];
    case 'y'
        outputArg = [0, -1i; 1i, 0];
    case 3
        outputArg = [1, 0; 0, -1];
    case 'z'
        outputArg = [1, 0; 0, -1];
    otherwise
        error('Undefined input argument for Pauli!')
end