import json
import os.path


def resolveFileInputJson(battmoDir, jsoninput):
    if isinstance(jsoninput, dict):
        if "isFile" in jsoninput:
            filename = jsoninput["filename"]
            fullfilename = os.path.join(battmoDir, filename)
            with open(fullfilename) as file:
                fileinput = json.load(file)
            jsoninput = fileinput
        for key in jsoninput:
            jsoninput[key] = resolveFileInputJson(battmoDir, jsoninput[key])
    return jsoninput


def loadJsonBattmo(battmoDir, filename):
    filename = os.path.join(battmoDir, filename)
    with open(filename) as file:
        jsoninput = json.load(file)
    jsoninput = resolveFileInputJson(battmoDir, jsoninput)
    return jsoninput
