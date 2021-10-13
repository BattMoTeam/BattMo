# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path

def batmoDir():
    return Path("/home/xavier/Matlab/Projects/project-batman/")

schema_folder = batmoDir() / Path('Liquid/CarbonateBased/')
schema_filename = schema_folder / 'binaryelectrolyte.schema.json'

base_uri = 'file://batmo/schemas/'

with open(schema_filename) as schema_file:
    schema = json.load(schema_file)

with open(instance_filename) as instance_file:
    instance = json.load(instance_file)

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
