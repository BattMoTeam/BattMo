import json
import inspect
from pathlib import Path
import os.path


def getBattMoDir():
    filename = inspect.getframeinfo(inspect.currentframe()).filename
    path = Path(os.path.dirname(os.path.abspath(filename)))
    path = path.parent.parent
    return path


def resolveFileInputJson(jsoninput):
    if type(jsoninput) is dict:
        if "isFile" in jsoninput:
            filename = jsoninput["filename"]
            fullfilename = getBattMoDir() / filename
            with open(fullfilename) as file:
                fileinput = json.load(file)
            jsoninput = fileinput
        for key in jsoninput:
            jsoninput[key] = resolveFileInputJson(jsoninput[key])
    return jsoninput


def loadJsonBatmo(filename):
    filename = getBattMoDir() / filename
    with open(filename) as file:
        jsoninput = json.load(file)
    jsoninput = resolveFileInputJson(jsoninput)
    return jsoninput
