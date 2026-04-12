#!/usr/bin/env python

# Copyright 2014  Steve Conklin 
# steve at conklinhouse dot com
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# You will need the (very nice) pySerial module, found here:
# http://pyserial.wiki.sourceforge.net/pySerial

import sys
#import os
#import os.path
#import string
#from array import *

if len(sys.argv) != 2:
    print 'Usage: %s filename' % sys.argv[0]
    sys.exit()

with open(sys.argv[1], "r") as nec:
    lns = nec.readlines()

print "Tag,Element,Segments, X1, Y1, Z1, X2, Y2, Z2, Radius"
for line in lns:
    if line.startswith("GW"):
        tag, element, numsegs, x1, y1, z1, x2, y2, z2, radius = line.split()
        print "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (tag, element, numsegs, x1, y1, z1, x2, y2, z2, radius)
