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

if len(sys.argv) != 3:
    print 'Usage: %s filename scalefactor' % sys.argv[0]
    sys.exit()

scalefactor = float(sys.argv[2])

with open(sys.argv[1], "r") as nec:
    lns = nec.readlines()

for line in lns:
    if line.startswith("GW"):
        tag, element, numsegs, x1, y1, z1, x2, y2, z2, radius = line.split()
        x1f = float(x1)
        y1f = float(y1)
        z1f = float(z1)
        x2f = float(x2)
        y2f = float(y2)
        z2f = float(z2)

        # what is the orientation?
        # If Z is zero, then elements are along Y
        # if Y is zero, then elements are along Z
        if (abs(z1f) < .001):
            # Along Y
            y1f = y1f * scalefactor
            y2f = y2f * scalefactor
        else:
            # Along Z
            z1f = z1f * scalefactor
            z2f = z2f * scalefactor
        print '{:5s} {:5s} {:2s} {: 3.5E} {: 3.5E} {: 3.5E} {: 3.5E} {: 3.5E} {: 3.5E} {:13s}'.format(tag, element, numsegs, x1f, y1f, z1f, x2f, y2f, z2f, radius)
        #print "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (tag, element, numsegs, x1, y1, z1, x2, y2, z2, radius)
    else:
        print line.strip()
