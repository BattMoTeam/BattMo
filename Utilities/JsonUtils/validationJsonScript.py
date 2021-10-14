# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson


def batmoDir():
    return Path("/home/xavier/Matlab/Projects/project-batman/")


schema_folder = batmoDir() / Path('JsonSchemas')
schema_filename = schema_folder / 'electrolyte.schema.json'

base_uri = 'file://batmo/schemas/'

with open(schema_filename) as schema_file:
    mainschema = json.load(schema_file)

# example
jsoninput = rjson.loadJsonBatmo('Battery/lithiumbattery.json')

resolver = jsonschema.RefResolver(base_uri=base_uri, referrer=schema)


def addJsonSchema(jsonschemaName):
    jsonchemaFilename = jsonschemaName + '.schema.json'
    schema_filename = schema_folder / jsonchemaFilename
    with open(schema_filename) as schema_file:
        refschema = json.load(schema_file)
    key = "file://batmo/schemas/" + jsonschemaName
    resolver.store[key] = refschema


schemaList = ["activematerial", "battery", "binaryelectrolyte",
              "currentcollector", "electrodeactivecomponent",
              "electrode", "electrolyte", "separator", "thermalmodel"]
for schema in schemaList:
    try:
        addJsonSchema(schema)
    except:
        print(schema)

v = jsonschema.Draft7Validator(mainschema, resolver=resolver)

v.is_valid(jsoninput)
