# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson

def batmoDir():
    return Path("/home/xavier/Matlab/Projects/project-batman/")

schema_folder = batmoDir() / Path('Materials/Liquid/CarbonateBased/')
schema_filename = schema_folder / 'binaryelectrolyte.schema.json'

base_uri = 'file://batmo/schemas/'

with open(schema_filename) as schema_file:
    schema = json.load(schema_file)

# example
jsoninput = rjson.loadJsonBatmo('Materials/Liquid/CarbonateBased/orgLiPF6.json')

resolver = jsonschema.RefResolver(base_uri=base_uri, referrer=schema)

schema_filename = schema_folder / 'electrolyte.schema.json'
with open(schema_filename) as schema_file:
    refschema = json.load(schema_file)
resolver.store["file://batmo/schemas/electrolyte"] = refschema

schema_filename = schema_folder / 'separator.schema.json'
with open(schema_filename) as schema_file:
    refschema = json.load(schema_file)
resolver.store["file://batmo/schemas/separator"] = refschema

v = jsonschema.Draft7Validator(schema, resolver=resolver)

v.is_valid(jsoninput)
