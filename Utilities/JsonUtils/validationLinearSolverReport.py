# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson
import os

schema_folder = rjson.getBattMoDir() / Path('Utilities') / Path('JsonSchemas')

# add "linesrsolver.schema" in resolver
base_uri = 'file://battmo/schemas/'
resolver = jsonschema.RefResolver(base_uri=base_uri, referrer={})
jsonSchemaFilename = Path('linearsolver.schema.json')
schema_filename = schema_folder / jsonSchemaFilename
with open(schema_filename) as schema_file:
    refschema = json.load(schema_file)

resolver.store["file://battmo/schemas/linearsolver"] = refschema

# load main schema
schema_filename = schema_folder / 'linearsolverreport.schema.json'
with open(schema_filename) as schema_file:
    mainschema = json.load(schema_file)

resolver.store["file://battmo/schemas/linearsolverreport"] = mainschema

v = jsonschema.Draft202012Validator(mainschema, resolver=resolver)

# jsonfilenames = ['linearsolver1.json',
                 # 'linearsolver2.json',
                 # 'linearsolver3.json']

jsonfilenames = ['linearSolverReportTest.json']

for jsonfilename in jsonfilenames:
    jsonfilename = schema_folder / 'Tests' / jsonfilename
    with open(jsonfilename) as jsonfile:
        jsoninput = json.load(jsonfile)
    v.validate(jsoninput)
    if v.is_valid(jsoninput):
        print('ok')
