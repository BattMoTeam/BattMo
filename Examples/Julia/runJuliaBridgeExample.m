jsonfilename = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', 'p2d_40.json');
jsonstruct   = parseBattmoJson(jsonfilename);

% The server should be start only once.
man = ServerManager();

output = man.runBatteryJson(jsonstruct)
