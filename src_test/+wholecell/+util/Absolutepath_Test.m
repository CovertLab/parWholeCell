classdef Absolutepath_Test < TestCase
    methods
        function this = Absolutepath_Test(name)
            this = this@TestCase(name);
        end
        
        function testAbsolutepath(~)
            tmp = pwd;
            if ispc
                tmp(1) = upper(tmp(1));
            end
            assertEqual(tmp(1:find(tmp == filesep, 1, 'last') - 1),  absolutepath('..'));
            assertEqual(pwd,  absolutepath('.'));
            assertEqual(fullfile(pwd, 'runTests.m'),  absolutepath('./runTests.m'));
        end
    end
end