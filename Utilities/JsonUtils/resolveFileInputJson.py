import json
from pathlib import Path

def batmoDir():
    return Path("/home/xavier/Matlab/Projects/project-batman/")

def resolveFileInputJson(jsoninput):
    if type(jsoninput) is dict:
        if "isFile" in jsoninput:
            filename = jsoninput["filename"]
            fullfilename = batmoDir() / filename
            with open(fullfilename) as file:
                fileinput = json.load(file)
            jsoninput = fileinput
        for key in jsoninput:
            jsoninput[key] = resolveFileInputJson(jsoninput[key])
    return jsoninput


def loadJsonBatmo(filename):
    filename = batmoDir() / filename
    with open(filename) as file:
        jsoninput = json.load(file)
    jsoninput = resolveFileInputJson(jsoninput)
    return jsoninput
