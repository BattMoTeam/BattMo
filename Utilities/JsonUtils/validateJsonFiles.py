# Call from MATLAB as py.validationJsonScript.validate(filename).
# Remember to restart the pyenv for changes to be seen by MATLAB. This
# is done by calling `terminate(pyenv)` or using the loadModule
# script.

import json
import os
import sys
from pathlib import Path
import jsonschema
from referencing import Registry, Resource
import resolveFileInputJson as rjson

VERBOSE = True


def addSchema(schemaFolder, registry, schemaFilename):
    """Add the jsonschema in the registry"""
    schemaFullFile = schemaFolder / schemaFilename
    if VERBOSE:
        print(schemaFullFile)
    with open(schemaFullFile) as fid:
        schema = json.load(fid)
    resource = Resource.from_contents(schema)
    baseUri = "file://./"
    uri = baseUri + schemaFilename
    registry = registry.with_resource(uri=uri, resource=resource)
    return registry


def validate(battmoDir, jsonfile, schemaFilename="Simulation.schema.json"):
    schemaFolder = battmoDir / Path("Utilities") / Path("JsonSchemas")
    if VERBOSE:
        print(f"{schemaFolder=}")

    # We collect the schema.json files and add them in the resolver
    registry = Registry()
    for _, _, filenames in os.walk(schemaFolder):
        for filename in filenames:
            if filename.endswith(".schema.json") and filename != schemaFilename:
                registry = addSchema(schemaFolder, registry, filename)

    # We validate the battery schema
    schemaFullFile = schemaFolder / schemaFilename
    with open(schemaFullFile) as fid:
        mainschema = json.load(fid)
    if VERBOSE:
        print("Validate main schema", mainschema)
    validator = jsonschema.Draft202012Validator(mainschema, registry=registry)

    # Validate the input jsonfile
    if VERBOSE:
        print("Validate input file", jsonfile)
    jsonstruct = rjson.loadJsonBattmo(battmoDir, jsonfile)
    validator.validate(jsonstruct)

    return True


if __name__ == "__main__":
    # Allow for command line call
    battmoDir = sys.argv[1]
    jsonfile = sys.argv[2]

    if len(sys.argv) == 3:
        validate(battmoDir, jsonfile)
    elif len(sys.argv) == 4:
        schemaFilename = sys.argv[3]
        validate(battmoDir, jsonfile, schemaFilename)
    else:
        raise RuntimeError(
            "Use either 2 or 3 arguments, right now it was", len(sys.argv)
        )
