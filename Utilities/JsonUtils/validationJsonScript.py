# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson
import os

schema_folder = rjson.getBattMoDir() / Path('Utilities') / Path('JsonSchemas')

base_uri = 'file://batmo/schemas/'

resolver = jsonschema.RefResolver(base_uri=base_uri, referrer={})


def addJsonSchema(jsonSchemaName, verbose=False):
    """Add the jsonschema in the resolver"""
    jsonSchemaFilename = jsonSchemaName + '.schema.json'
    schema_filename = schema_folder / jsonSchemaFilename
    if verbose:
        print(schema_filename)
    with open(schema_filename) as schema_file:
        refschema = json.load(schema_file)
    key = "file://batmo/schemas/" + jsonSchemaName
    resolver.store[key] = refschema


# We collect the schema.json files and add them in the resolver
for (dirpath, dirnames, filenames) in os.walk(schema_folder):
    for filename in filenames:
        if filename.endswith('.schema.json'):
            jsonSchemaName = filename.replace('.schema.json', '')
            addJsonSchema(jsonSchemaName, verbose=True)

# We validate the battery schema
schema_filename = schema_folder / 'Battery.schema.json'
with open(schema_filename) as schema_file:
    mainschema = json.load(schema_file)

v = jsonschema.Draft202012Validator(mainschema, resolver=resolver)

jsonfiles = ['ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json',
             'ParameterData/ParameterSets/Xu2015/lfp.json',
             'ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json']

for jsonfile in jsonfiles:
    print(jsonfile)
    jsoninput = rjson.loadJsonBatmo(jsonfile)
    v.validate(jsoninput)
    if v.is_valid(jsoninput):
        print('ok')
