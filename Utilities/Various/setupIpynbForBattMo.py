import sys
import re

def fixBattMoCells(inputstr):
    """
    Search each cell and make sure the cell type is set correctly in the cells that contain a battmo link

    Parameters:
    inputstr (str): string that has been obtained by parsing a ipynb notebook

    Returns:
    str : input string where the cells in markdown have been changed to raw format, which can be processes by nbsphinx
    """

    dofix = inputstr.find(':battmo', 0)
    while dofix > 0:
        [inputstr, pos] = switchCellToRaw(inputstr, dofix)
        dofix = inputstr.find(':battmo', pos)
    return inputstr

def switchCellToRaw(inputstr, fix_position):

    """Starts at the fix_position in the string inputstr. From there, finds the containing cell and changes its type and
    metadata mimetype. Returns the modified inputstr and the position of the end of the cell.
    """
    
    cell_start_pos = inputstr.rfind('"cell_type"', 0, fix_position)
    pos = inputstr.find('"cell_type": "markdown"', cell_start_pos)

    if pos < fix_position:
        inputstr = inputstr[:cell_start_pos] + inputstr[cell_start_pos:].replace('markdown', 'raw', 1)
        # inputstr = inputstr.replace('markdown', 'raw')
        match = re.search('"metadata": {', inputstr[cell_start_pos:-1])
        pos = cell_start_pos + match.end()
        s = '\n"raw_mimetype": "text/restructuredtext"\n'
        inputstr = inputstr[:pos] + s + inputstr[pos:]

    cell_start_pos = inputstr.rfind('{', 0, cell_start_pos)

    cell_end_pos = getMatchingBracketPosition(inputstr, cell_start_pos)

    return inputstr, cell_end_pos

def getMatchingBracketPosition(inputstr, startpos):

    """ Utility function to find matching parenthesis which is used to find the end position of a cell"""
    
    inputstr = inputstr[startpos:]
    count = 0
    for i, s in enumerate(inputstr):
        if s == '}':
            count -= 1
        elif s == '{':
            count += 1
        if count == 0:
            return startpos + i + 1
        
        

if __name__ == "__main__":

    if len(sys.argv) == 1:
        raise ValueError("No input file provided. Please provide the path to the input file as an argument.")
    else:
        input_filename = sys.argv[1]
        if len(sys.argv) == 3:
            output_filename = sys.argv[2]
        else:
            output_filename = input_filename
            
    with open(input_filename, 'r') as file:
        inputstr = file.read()

    inputstr = fixBattMoCells(inputstr)

    with open(output_filename, 'w') as file:
        file.write(inputstr)
