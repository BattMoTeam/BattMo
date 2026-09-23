import json
import os.path


def resolveFileInputJson(battmoDir, jsoninput):
    if isinstance(jsoninput, list):
        return [resolveFileInputJson(battmoDir, item) for item in jsoninput]

    if isinstance(jsoninput, dict):
        if "isFile" in jsoninput:
            filename = jsoninput["filename"]
            fullfilename = os.path.join(battmoDir, filename)
            with open(fullfilename) as file:
                fileinput = json.load(file)
            return resolveFileInputJson(battmoDir, fileinput)

        return {
            key: resolveFileInputJson(battmoDir, value)
            for key, value in jsoninput.items()
        }

    return jsoninput


def mergeJson(primary, secondary):
    """Recursively merge dictionaries, keeping values from primary on conflicts."""
    result = dict(primary)
    for key, value in secondary.items():
        if key not in result:
            result[key] = value
        elif isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = mergeJson(result[key], value)
    return result


def resolveMergeInputJson(jsoninput):
    if isinstance(jsoninput, list):
        return [resolveMergeInputJson(item) for item in jsoninput]

    if not isinstance(jsoninput, dict):
        return jsoninput

    if "MergeInputs" in jsoninput:
        merge_input = jsoninput["MergeInputs"]
        inputs = [resolveMergeInputJson(item) for item in merge_input["Inputs"]]
        if merge_input.get("merge_order", "ascending") == "ascending":
            inputs.reverse()

        result = {}
        for item in inputs:
            result = mergeJson(result, item)
        return result

    return {
        key: resolveMergeInputJson(value)
        for key, value in jsoninput.items()
    }


def loadJsonBattmo(battmoDir, filename):
    filename = os.path.join(battmoDir, filename)
    with open(filename) as file:
        jsoninput = json.load(file)
    jsoninput = resolveFileInputJson(battmoDir, jsoninput)
    jsoninput = resolveMergeInputJson(jsoninput)
    return jsoninput
