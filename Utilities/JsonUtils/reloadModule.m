function reloadModule(modulename)

    if nargin == 0
        modulename = 'validationJsonScript';
    end

    mod = py.importlib.import_module(modulename);
    py.importlib.reload(mod);

end
