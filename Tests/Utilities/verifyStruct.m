function verifyStruct(testCase, state, testname)

    filename = sprintf('%s.json', testname);
    refstate = parseBattmoJson(fullfile('Tests', 'TestExamples', 'ReferenceData', filename));

    import matlab.unittest.constraints.IsEqualTo
    import matlab.unittest.constraints.RelativeTolerance

    % Compare each data point in state with refstate
    reltol = 1e-6;
    testCase.assertThat(state, IsEqualTo(refstate, ...
                                         'Within', RelativeTolerance(reltol)));

end
