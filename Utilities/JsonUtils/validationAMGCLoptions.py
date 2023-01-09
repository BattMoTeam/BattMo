# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson
import os

schema_folder = rjson.getBattMoDir() / Path('Utilities') / Path('JsonSchemas')

schema_filename = schema_folder / 'amgclOptions.schema.json'
with open(schema_filename) as schema_file:
    mainschema = json.load(schema_file)

v = jsonschema.Draft202012Validator(mainschema)

jsonfilenames = ['amgcl_options.json']

for jsonfilename in jsonfilenames:
    jsonfilename = schema_folder / 'Tests' / jsonfilename
    with open(jsonfilename) as jsonfile:
        jsoninput = json.load(jsonfile)
    v.validate(jsoninput)
    if v.is_valid(jsoninput):
        print('ok')
