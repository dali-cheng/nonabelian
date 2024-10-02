function outputArg = ternary(myCond, trueResult, falseResult)
%TERNARY Ternary operator in MATLAB
%   Returns the value (myCond) ? (trueResult) : (falseResult).
if myCond
    outputArg = trueResult;
else
    outputArg = falseResult;
end
end