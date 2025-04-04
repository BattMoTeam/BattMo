# Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
# 
# This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
# 
# MRST is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# MRST is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with MRST.  If not, see <http://www.gnu.org/licenses/>.

import os
import xml.etree.ElementTree
import re
import sys
from os.path import join, splitext, exists

publishedExampleDir = os.path.realpath(__file__)
publishedExampleDir = os.path.dirname(publishedExampleDir)
publishedExampleDir = join(publishedExampleDir, '..', 'publishedExamples')

def parse_publish_to_xml(exfile):
    tmp = os.path.split(exfile)
    example_filename = tmp[-1].replace('.xml','.m')
    f = open(exfile, 'r')
    example = f.read()
    f.close()
    example = example.replace("<tt>", "`")
    example = example.replace("</tt>", "`")
    # example = re.sub("<a [^>]*>([^<]*)</a>", "", example)
    example = re.sub(r'<a [^>]*>([^<]*)</a>', r'XMLLT\1XMLGT', example)
    #example = re.sub('\.\.\.\s+', '...' + os.linesep, example)
    #example = example.replace("...", "...\n")
    e = xml.etree.ElementTree.fromstring(example)
    # e = xml.etree.ElementTree.parse(exfile).getroot()
    # outfile  = open("D:/jobb/bitbucket/mrst-documentation/flowSolverTutorial1.txt", "w")
    parsedText = ''
    rstDescriptionFile = tmp[-1].replace('.xml', 'Preamble.rst')
    if exists(join(publishedExampleDir, rstDescriptionFile)):
        includeRstDescriptionFile = ".. include:: " + rstDescriptionFile + "\n\n"
    else:
        includeRstDescriptionFile = "\n"

    # include link
    examplename = tmp[-1].replace('.xml', '')

    def addIndent(txt):
        return "  ".join(("\n"+txt.lstrip()).splitlines(True)) + "\n"

    headcount = 0
    for el in e:
        if el.tag == "cell":
            # A cell which corresponds to a matlab section in cell mode
            for paragraphs in el:
                ptag = paragraphs.tag;
                if ptag == "steptitle":
                    # Heading of subsection
                    if headcount == 0:
                        sect =  ".. _" + examplename + ":\n\n"
                    else:
                        sect = ""
                    tit = "".join(paragraphs.itertext());
                    tit_length = len(tit);
                    if headcount == 0:
                        # Make heading
                        ch = '='
                        sect = sect + '\n' + ch * tit_length + '\n' + tit + '\n' + ch * tit_length + '\n'
                    else:
                        ch = '='
                        sect = sect + '\n' + tit + '\n' + ch * tit_length + '\n'                    
                    if headcount == 0:
                        sect = sect + '*Generated from ' + example_filename + '*\n\n'
                        sect = sect + includeRstDescriptionFile
                        headcount += 1;
                    parsedText += sect
                elif ptag == "text":
                    # Simply text
                    for block in paragraphs:
                        if block.tag == "p":
                            # We found text
                            if block.text.strip():
                                #txt = block.text
                                #parsedText += txt
                                #parsedText += "\n"
                                parsedText += "".join(block.itertext());
                            else:
                                # This paragraph is an equation
                                for subblock in block:
                                    if subblock.tag == "equation":
                                        parsedText += "\n.. math::\n"
                                        eq = subblock.get("text")
                                        if eq:
                                            eq = eq.strip("$")
                                            eq = addIndent(eq)
                                            parsedText += eq
                        parsedText += "\n"
                elif ptag == "mcode" or ptag == "mcode-xmlized":
                    # Code block
                    for block in paragraphs:
                        #  print block
                        txt = "".join(block.itertext());
                        txt = re.sub(r'\.\.\.\s', '...' + os.linesep, txt)

                        crindex = txt.lower().find("copyright")
                        if crindex > -1:
                            # Strip copyright plus any comments etc before that.
                            txt = txt[0:crindex]
                            lastalpha = [m.end(0) for m in re.finditer(r'\w|\)|\;', txt)];
                            if lastalpha:
                               txt = txt[0:lastalpha[-1]]
                        parsedText += "\n.. code-block:: matlab\n"
                        txt = addIndent(txt);
                        parsedText += txt
                        parsedText += "\n"
                elif ptag == "mcodeparsedText":
                    # Another code block?? Slightly different format
                    parsedText += "\n.. code-block:: none\n"
                    txt = addIndent(paragraphs.text);
                    parsedText += txt
                    parsedText += "\n"
                elif ptag == "img":
                    # Image. Use provided relative path.
                    img = paragraphs.get("src")
                    parsedText += '.. figure:: ' + img + "\n"
                    parsedText += "  :figwidth: 100%\n"
                    parsedText += "\n"
        elif el.tag == "originalCode":
            sourcecode = ":orphan:\n\n"
            sourcecode += ".. _" + examplename + "_source:\n\n"
            title = "Source code for " + examplename
            sourcecode += title + "\n"
            sourcecode += len(title)*"-" + "\n\n"
            sourcecode += ".. code:: matlab\n\n"
            sourcecode += addIndent(el.text)

    parsedText = parsedText.replace("XMLLT", "<")
    parsedText = parsedText.replace("XMLGT", ">")

    parsedText += "\n\n"
    parsedText += "complete source code can be found :ref:`here<" + examplename + "_source>`\n"
    output = dict()
    output['parsed'] = parsedText
    output['source'] = sourcecode
    return output



if len(sys.argv) < 2:
    allfiles = os.walk(publishedExampleDir)
    allfiles = [f[-1] for f in allfiles]
    xmlfiles = [f for f in allfiles[0] if splitext(f)[-1] == '.xml']
else:
    xmlfiles = [join(publishedExampleDir, sys.argv[1] + '.xml')]

for xmlfile in xmlfiles:
    xmlpath = join(publishedExampleDir, xmlfile)
    examplename = splitext(xmlfile)[0]
    output = parse_publish_to_xml(xmlpath)
    outfile = open(join(publishedExampleDir, examplename + '.rst'), 'w')
    outfile.write(output['parsed'])
    outfile.close()
    outfile = open(join(publishedExampleDir, examplename + '_source.rst'), 'w')
    outfile.write(output['source'])
    outfile.close()
    print("parsed {}".format(xmlfile))
