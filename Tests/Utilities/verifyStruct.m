function verifyStruct(testCase, state, testname)

    filename = sprintf('%s.json', testname);
    refstate = parseBattmoJson(fullfile('Tests', 'TestExamples', 'ReferenceData', filename));

    import matlab.unittest.constraints.IsEqualTo

    % FIXME Add tolerance
    testCase.verifyThat(state, IsEqualTo(refstate));

end
