function o = callFunction(a,b)

    o = 0;
    
    % Call the terminate function
    o = coder.ceval('sum_operation',a,b);
    
 
end