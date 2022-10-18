# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson


schema_folder = rjson.getBattMoDir() / Path('Utilities') / Path('JsonSchemas')

base_uri = 'file://batmo/schemas/'

resolver = jsonschema.RefResolver(base_uri=base_uri, referrer={})


def addJsonSchema(jsonschemaName):
    jsonchemaFilename = jsonschemaName + '.schema.json'
    schema_filename = schema_folder / jsonchemaFilename
    with open(schema_filename) as schema_file:
        refschema = json.load(schema_file)
    key = "file://batmo/schemas/" + jsonschemaName
    resolver.store[key] = refschema


schemaList = ["activematerial",
              "battery",
              "binaryelectrolyte",
              "interface",
              "currentcollector",
              "soliddiffusion",
              "electrode",
              "electrolyte",
              "separator",
              "thermalmodel"]

for schema in schemaList:
    addJsonSchema(schema)

# We validate the battery schema
schema_filename = schema_folder / 'battery.schema.json'
with open(schema_filename) as schema_file:
    mainschema = json.load(schema_file)

v = jsonschema.Draft7Validator(mainschema, resolver=resolver)

jsoninput = rjson.loadJsonBatmo('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json')

if v.is_valid(jsoninput):
    print('ok')
else:
    print('not ok')
