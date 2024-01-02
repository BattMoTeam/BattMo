# Call from MATLAB as py.validationJsonScript(filename). Remember to
# restart the pyenv for changes to be seen by MATLAB. This is done by
# calling `terminate(pyenv)` or using the reloadModule script.

import json
import jsonschema
from pathlib import Path
import resolveFileInputJson as rjson
import os
import sys

schema_folder = rjson.getBattMoDir() / Path("Utilities") / Path("JsonSchemas")
base_uri = "file://./"
resolver = jsonschema.RefResolver(base_uri=base_uri, referrer={})
verbose = False


def addJsonSchema(jsonSchemaName):
    """Add the jsonschema in the resolver"""
    jsonSchemaFilename = jsonSchemaName + ".schema.json"
    schema_filename = schema_folder / jsonSchemaFilename
    if verbose:
        print(schema_filename)
    with open(schema_filename) as schema_file:
        refschema = json.load(schema_file)
    key = base_uri + jsonSchemaFilename
    resolver.store[key] = refschema


def validate(jsonfile):
    # We collect the schema.json files and add them in the resolver
    for dirpath, dirnames, filenames in os.walk(schema_folder):
        for filename in filenames:
            if filename.endswith(".schema.json"):
                jsonSchemaName = filename.replace(".schema.json", "")
                addJsonSchema(jsonSchemaName)

    # We validate the battery schema
    schema_filename = schema_folder / "Simulation.schema.json"
    with open(schema_filename) as schema_file:
        mainschema = json.load(schema_file)
    if verbose:
        print("Validate main schema", mainschema)
    v = jsonschema.Draft202012Validator(mainschema, resolver=resolver)

    # Validate the input jsonfile
    if verbose:
        print("Validate input file", jsonfile)
    jsonstruct = rjson.loadJsonBattmo(jsonfile)
    v.validate(jsonstruct)

    return True


if __name__ == "__main__":
    # Allow for command line call
    validate(sys.argv[1])
