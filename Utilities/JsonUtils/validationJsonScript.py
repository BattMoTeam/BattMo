# Example usage for add_local_schemas_to

import json
import jsonschema
from pathlib import Path
from urllib.parse import urljoin

schema_folder = Path('/home/xavier/Matlab/Projects/project-batman/Materials/Liquid/CarbonateBased/')
schema_filename = schema_folder / 'binaryelectrolyte.schema.json'
instance_folder = schema_folder
instance_filename = instance_folder / 'orgLiPF6.json'
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

v = jsonschema.Draft7Validator(schema, resolver=resolver)

v.is_valid(instance)
