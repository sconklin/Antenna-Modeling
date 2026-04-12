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
import serial
from telnetlib import *
import json
import time

class foxdelta():

    def __init__(self):
        self.verbose = True
        self.ser = None
        return

    def open(self, cport='/dev/ttyUSB0', baud=9600):
        self.ser = serial.Serial(port=cport, baudrate=baud, parity='N', stopbits=1, timeout=3, xonxoff=0, rtscts=0, dsrdtr=0)
        if self.ser == None:
            print 'Unable to open serial device %s' % cport
            raise IOError
        return

    def close(self):
        if self.ser:
            self.ser.close()
        return

    def readline(self):
        inch = '0'
        line = ''
        while (inch != '\r'):
            inch = self.ser.read()
            line = line + inch
        return line

    def write(self, outbuf):
        self.ser.write(outbuf)
        return

    def getVersions(self):
        self.write(':99\r')
        vers = self.readline()
        #":99" ch5(AUTHOR_ID) ch10(PROJECT_ID) ch5(FW_ID) CR
        return {'author':vers[3:8], 'project':vers[9:19], 'firmware':vers[20:25]}

    def makeSingleMeasurement(self, frequency):
        #       msgOUT = ":24" dec4(vIn0) dec4(RLdB*10) dec4(SWR*100) dec9(Frequency) S(+/-)Dec4(dB0*100) CR
        self.write(':24%09d\r' % frequency)
        meas = self.readline()
        #print meas
        results =  {'vIn0': int(meas[3:7]), 'RLdB': int(meas[7:11])/ 10.0, \
                    'SWR': int(meas[11:15])/100.0, 'Freq': int(meas[15:24]), \
                    'dB0': int(meas[24:28])/100.0}
        #print "%d : %f" % (results['Freq'], results['SWR'])
        return results

    def makeMultipleMeasurements(self, startFreq, endFreq, freqStep):
        results = []
        for freq in range(startFreq, endFreq, freqStep):
            results.append(self.makeSingleMeasurement(freq))
        return results

class rotator():

    def __init__(self):
        self.verbose = True
        self.tel = None
        return

    def open(self, ip='172.31.0.15', port=4533):
        self.tel = Telnet(ip, port)
        if self.tel == None:
            print 'Unable to open telnet connection to %s:%d' % (ip, port)
            raise IOError
        return

    def close(self):
        if self.tel:
            self.tel.close()
        return

    def dumplines(self):
        # read until we time out
        while True:
            line = self.readline()
            if line == '':
                return
            #else:
            #    print "dumping line %s" % line

    def readline(self):
        line = self.tel.read_until('\n', 3)
        return line

    def write(self, outbuf):
        self.tel.write(outbuf)
        return

    def setAzimuth(self, azimuth):
        #print "setAzimuth %d" % azimuth
        self.write('P %d 0\n' % azimuth)
        self.dumplines()
        return

    def getAzimuth(self):
        #print "getAzimuth"
        self.write('p\n')
        az = self.readline()
        el = self.readline()
        return float(az.strip())

if len(sys.argv) != 2:
    print 'Usage: %s serialdevice' % sys.argv[0]
    sys.exit()

freqStep = 1000
freqList = [
#    {'band': '20m', 'start':14000000, 'end':14350000},
#    {'band': '17m', 'start':18068000, 'end':18168000},
    {'band': '17m', 'start':18070000, 'end':18167000},
#    {'band': '15m', 'start':21000000, 'end':21700000},
#    {'band': '12m', 'start':24890000, 'end':24990000},
#    {'band': '10m', 'start':28000000, 'end':29700000},
#    {'band': '6m', 'start':50000000, 'end':54000000}
]

# Open the SWR measuring device
fd = foxdelta()
fd.open(cport=sys.argv[1], baud=38400)

# Open the rotator control
rc = rotator()
rc.open(ip="172.31.0.15", port=4533)

lastAz = 0
firstTime = True

#Where do we keep the data?
alldata = {}
for f in freqList:
    alldata[f['band']] = []

#for az in range(0, 20, 5):
for az in range(0, 360, 5):
    print "Setting Azimuth = %d" % az
    # Check for the full rotation case and allow time
    rc.setAzimuth(az)
    if (((lastAz <=180) and (az > 180)) or firstTime):
        firstTime = False
        time.sleep(60)
    else:
        time.sleep(5)
    lastAz = az
    # sometimes it overshoots
    rc.setAzimuth(az)
    time.sleep(5)

    actualAz = rc.getAzimuth()
    print "Read Actual Azimuth = %s" % actualAz
    
    for f in freqList:
        print "Testing %s band" % f['band']
        result = fd.makeMultipleMeasurements(f['start'], f['end']+1, freqStep)
        out = {}
        out['band'] = f['band']
        out['azimuth'] = actualAz
        out['data'] = result
        alldata[f['band']].append(out)
        
for f in freqList:
    fname = "%s.json" % f['band']
    with open(fname, "w") as fo:
        fo.write(json.dumps(alldata[f['band']], indent=4))

fd.close()
rc.close()
