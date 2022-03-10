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
from os.path import *

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
    output = ''
    rstDescriptionFile = tmp[-1].replace('.xml', 'Preamble.rst')
    includeRstDescriptionFile = ".. include:: " + rstDescriptionFile + "\n\n"

    def addIdent(txt):
        return "  ".join(("\n"+txt.lstrip()).splitlines(True)) + "\n"

    headcount = 0
    for el in e:
        if el.tag == "cell":
            # A cell which corresponds to a matlab section in cell mode
            for paragraphs in el:
                ptag = paragraphs.tag;
                if ptag == "steptitle":
                    # Heading of subsection
                    sect = "".join(paragraphs.itertext());
                    tit_length = len(sect);
                    if headcount == 0:
                        # Make heading
                        ch = '-';

                    else:
                        ch = '^'
                    sect = '\n' + sect + '\n' + ch * tit_length + '\n'
                    if headcount == 0:
                        sect = sect + '*Generated from ' + example_filename + '*\n\n'
                        sect = sect + includeRstDescriptionFile
                    headcount += 1;
                    output += sect
                elif ptag == "text":
                    # Simply text
                    for block in paragraphs:
                        if block.tag == "p":
                            # We found text
                            if block.text.strip():
                                #txt = block.text
                                #output += txt
                                #output += "\n"
                                output += "".join(block.itertext());
                            else:
                                # This paragraph is an equation
                                for subblock in block:
                                    if subblock.tag == "equation":
                                        output += "\n.. math::\n"
                                        eq = subblock.get("text")
                                        if eq:
                                            eq = eq.strip("$")
                                            eq = addIdent(eq)
                                            output += eq
                        output += "\n"
                elif ptag == "mcode" or ptag == "mcode-xmlized":
                    # Code block
                    for block in paragraphs:
                        #  print block
                        txt = "".join(block.itertext());
                        txt = re.sub('\.\.\.\s', '...' + os.linesep, txt)

                        crindex = txt.lower().find("copyright")
                        if crindex > -1:
                            # Strip copyright plus any comments etc before that.
                            txt = txt[0:crindex]
                            lastalpha = [m.end(0) for m in re.finditer('\w|\)|\;', txt)];
                            if lastalpha:
                               txt = txt[0:lastalpha[-1]]
                        output += "\n.. code-block:: matlab\n"
                        txt = addIdent(txt);
                        output += txt
                        output += "\n"
                elif ptag == "mcodeoutput":
                    # Another code block?? Slightly different format
                    output += "\n.. code-block:: none\n"
                    txt = addIdent(paragraphs.text);
                    output += txt
                    output += "\n"
                elif ptag == "img":
                    # Image. Use provided relative path.
                    img = paragraphs.get("src")
                    output += '.. figure:: ' + img + "\n"
                    output += "  :figwidth: 100%\n"
                    output += "\n"
    output = output.replace("XMLLT", "<")
    output = output.replace("XMLGT", ">")
    return output

publishedExampleDir = os.path.realpath(__file__)
publishedExampleDir = os.path.dirname(publishedExampleDir)
publishedExampleDir = join(publishedExampleDir, '..', 'publishedExamples')

allfiles = os.walk(publishedExampleDir)
allfiles = [f[-1] for f in allfiles]
xmlfiles = [f for f in allfiles[0] if splitext(f)[-1] == '.xml']

for xmlfile in xmlfiles:
    xmlpath = join(publishedExampleDir, xmlfile)
    examplename = splitext(xmlfile)[0]
    exampletxt = parse_publish_to_xml(xmlpath)
    outfile = open(join(publishedExampleDir, examplename + '.rst'), 'w')
    outfile.write(exampletxt)
    outfile.close()
