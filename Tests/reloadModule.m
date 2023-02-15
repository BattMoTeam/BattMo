function reloadModule(modulename)

    clear classes;

    if nargin == 0
        modulename = 'validationJsonScript';
    end

    mod = py.importlib.import_module(modulename);
    py.importlib.reload(mod);

end
