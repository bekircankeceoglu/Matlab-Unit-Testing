classdef TestFunction < matlab.unittest.TestCase


    methods ( Test )
        
        function lessthanTest( testCase )

            o = callFunction(3,5);


            testCase.verifyLessThanOrEqual( o, 10);
            
        end
        
        function equalTest( testCase )

            o = callFunction(3,5);


            testCase.verifyEqual( o, 8);
            
        end
        
    end
end