#!/usr/bin/env python
# -*- coding: utf-8 -*-


__author__ = "Raony Guimarães"
__copyright__ = "Copyright 2015 Group By Extensions"
__credits__ = ["Raony Guimarães"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Raony Guimarães"
__email__ = "raonyguimaraes@gmail.com"
__status__ = "Development"

#group_by_extensions.py

from optparse import OptionParser
import os, shutil

parser = OptionParser()
usage = "usage: %prog [options] -i directory"
parser = OptionParser(usage=usage)

parser.add_option("-i", dest="input_dir",
                  help="Directory to group files by", metavar="DIR")

    
(options, args) = parser.parse_args()

input_dir=options.input_dir
output_dir="grouped_files"

for (dir, _, files) in os.walk(input_dir):
    for f in files:
        path = os.path.join(dir, f)
        # print path
        fileName, fileExtension = os.path.splitext(path)
        basename = os.path.basename(path)
        # print fileName, fileExtension
        fileExtension = fileExtension.replace('.', '')
        dest = '%s/%s/%s' % (output_dir, fileExtension, basename)
        directory = '%s/%s/' % (output_dir, fileExtension)
        if not os.path.exists(directory):
            os.makedirs(directory)
        shutil.copyfile(path, dest)
