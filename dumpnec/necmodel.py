#!/usr/bin/python
import sys

# Classes for each card type
class Taper(object):
    """
    A definition of a tapered wire
    """
    def __init__(self, line):
        self.input = line
        parts = line.split()
        self.ratio = float(parts[3])
        self.radiusFirst = float(parts[4])
        self.radiusLast = float(parts[5])
        return

    def dump():
        print "GC Wire Taper"
        print '    Ratio to prior segment:   %f' % self.ratio
        print '    First Segment Radius:      %f m (%.2f in)' % (self.radiusFirst, inches(self.radiusFirst))
        print '    Last  Segment Radius:      %f m (%.2f in)' % (self.radiusLast, inches(self.radiusLast))

class Wire(object):
    """
    Wire definition for an NEC2 model
    """
    def __init__(self, line):

        self.input = line
        parts = line.split()
        self.tag = int(parts[1])
        self.segments  = int(parts[2])
        self.startX = float(parts[3])
        self.startY = float(parts[4])
        self.startZ = float(parts[5])
        self.endX = float(parts[6])
        self.endY = float(parts[7])
        self.endZ = float(parts[8])
        self.wireRadius = float(parts[9])
        self.tobj = None
        return

    def dump(self):
        print 'GW Wire Specification:'
        print '    Tag:           %d' % self.tag
        print '    Segments:      %d' % self.segments
        print '    Start x:       %f m (%.2f in)' % (self.startX, inches(self.startX))
        print '    Start y:       %f m (%.2f in)' % (self.startY, inches(self.startY))
        print '    Start z:       %f m (%.2f in)' % (self.startZ, inches(self.startZ))
        print '    End   x:       %f m (%.2f in)' % (self.endX, inches(self.endX))
        print '    End   y:       %f m (%.2f in)' % (self.endY, inches(self.endY))
        print '    End   z:       %f m (%.2f in)' % (self.endZ, inches(self.endZ))
        print '    Wire Radius:   %f m (%.2f in)' % (self.wireRadius, inches(wireRadius))
        return

    def addTaper(self, tobj):
        # add a taper definition, detect two in a row (error)
        # TODO test for floating point zero
        if self.wireRadius != 0.0:
            raise ValueError('Taper definition added to a wire which had a radius defined')
        self.tobj = tobj
        return

class Arc(object):
    """
    Wire Arc definitions for an NEC2 model
    """
    # -spc- TODO can arcs have tapers?
    def __init__(self, line):

        self.input = line
        parts = line.split()
        self.tag = int(parts[1])
        self.segments = int(parts[2])
        self.radius = float(parts[3])
        self.angleStart = float(parts[4])
        self.angleEnd = float(parts[5])
        self.wireRadius = float(parts[6])
        return

    def dump(self):
        print 'GA Wire arc specification:'
        print '    Tag:          %d' % self.tag
        print '    Segments:     %d' % self.segments
        print '    Arc Radius:   %f m (%.2f in)' % (self.radius, inches(self.radius))
        print '    Angle start:  %f degrees' % self.angleStart
        print '    Angle end:    %f degrees' % self.AngleEnd
        print '    Wire Radius:  %f' % self.wireRadius
        return

class Helix(object):
    """
    Helix definition for an NEC2 model
    """
    def __init__(self, line):

        self.input = line
        parts = line.split()
        self.tag = int(parts[1])
        self.segments = int(parts[2])
        self.spacing = float(parts[3])
        self.length = float(parts[4])
        self.radiusXZ0 = float(parts[5])
        self.radiusXZhl = float(parts[6])
        self.radiusYZ0 = float(parts[7])
        self.radiusYZhl = float(parts[8])
        self.radius = float(parts[9])
        self.tobj = None
        return

    dump():
        print 'GH Helix/Spiral specification:'
        if self.length == 0:
            print '    (Spiral)'
        elif self.length > 0:
            print '    (Right Handed Helix)'
        else:
            print '    (Left Handed Helix)'
        print '    Segment Tag:         %d' % self.tag
        print '    Segments:            %d' % self.segments
        print '    Turn Spacing:        %f m (%.2f in)' % (self.spacing, inches(self.spacing))
        print '    Helix Length HL:     %f m (%.2f in)' % (self.length, inches(self.length))
        print '    Radius in x @ z=0:   %f m (%.2f in)' % (self.radiusXZ0, inches(self.radiusXZ0))
        print '    Radius in y @ z=0:   %f m (%.2f in)' % (self.radiusYZ0, inches(self.radiusYZ0))
        print '    Radius in x @ z=HL:  %f m (%.2f in)' % (self.radiusXZhl, inches(self.radiusXZhl))
        print '    Radius in y @ z=HL:  %f m (%.2f in)' % (self.radiusYZhl, inches(self.radiusYZhl))
        print '    Wire Radius:         %f' % self.radius
        return

    def addTaper(self, tobj):
        # -spc- TODO Can a helix have a taper????
        # add a taper definition, detect two in a row (error)
        # TODO test for floating point zero
        if self.wireRadius != 0.0:
            raise ValueError('Taper definition added to a wire which had a radius defined')
        self.tobj = tobj
        return

class Cylinder(object):
    """
    Cylinder definition for an NEC2 model

    To reproduce a structure by rotating about the Z-axis to form a
    complete cylindrical array, and to set flags so that symmetry is
    utilized in the solution.
    """
    # The tag increment (I1) is used to avoid duplication of tag numbers in
    # the reproduced structures. In forming a new structure for the array,
    # all valid tags on the previous copy or original structure are
    # incremented by (I1). Tags equal to zero are not incremented.
    # * The GR card should never be used when there are segments on the Z-axis
    # or crossing the Z-axis since overlapping segments would result.

    def __init__(self, line):
        self.line  = line
        parts = line.split()
        self.tagIncrement = int(parts[1])
        self.numOccurrances = int(parts[2])
        # TODO do something here????
        return

    def dump(self):
            print 'GR Generate cylindrical structure:'
            print '    Tag Increment:                     %d' % self.tagIncrement
            print '    Number of times structure occurs:  %d' % self.numOccurrances
        return

class Transform(object):
    """
    Geometry transformation for an NEC2 model
    """
    def __init__(self, line):
        self.line  = line
        parts = line.split()
        self.tagIncrement =  int(parts[1])
        self.newStructs = int(parts[2])
        self.rotationX = float(parts[3])
        self.rotationY = float(parts[4])
        self.rotationZ = float(parts[5])
        self.translationX = float(parts[6])
        self.translationY = float(parts[7])
        self.translationZ = float(parts[8])
        # Initial Translation Segment, round to an integer
        self.its = float(parts[9])
        return

    def dump(self):
        print 'GM Coordinate transformation:'
        print '    Tag Increment:                %d' % self.tagIncrement
        print '    New Structs                   %d' % self.newStructs
        print '    Rotation about x:             %f degrees' % self.rotationX
        print '    Rotation about y:             %f degrees' % self.rotationY
        print '    Rotation about z:             %f degrees' % self.rotationZ
        print '    Translation in x:             %f m (%.2f in)' % (self.translationX, inches(self.translationX))
        print '    Translation in y:             %f m (%.2f in)' % (self.translationY, inches(self.translationY))
        print '    Translation in z:             %f m (%.2f in)' % (self.translationZ, inches(self.translationZ))
        print '    Initial Translation Segment:  %d' % int(round(self.its))
        return

class Scale(object):
    """
    Geometry scaling for an NEC2 model
    """
    def __init__(self, line):
        self.line  = line
        parts = line.split()
        self.factor = float(parts[3])
        return

    def dump(self):
        print 'GS Scale Structure Domensions:'
        print '    Scale Factor:                 %f' % self.factor
        return

class GeometryEnd(object):
    """
    The end of the geometry definition
    """
    def __init__(self, line):
        self.input = line
        parts = line.split()
        self.groundPlane = int(parts[1])
        return

    def dump():
        print '    Ground Plane Flag:   %d' % self.groundPlane, 
        if self.groundPlane == 0:
            print ' (No Ground Plane Present)'
        elif self.groundPlane == 1:
            print ' (Ground Plane - Current Expansion Modified)'
        elif self.groundPlane == -1:
            print ' (Ground Plane - Current Expansion Not Modified)'
        else:
            print ' (Unknown Value)'
        return

# Classes for each subsection
class Geometry(object):
    """
    Antenna geometry information
    """
    def __init(self):
        self.lastWireAdded = None
        self.wiresInOrder = []
        self.wires = {}
        self.transforms = []
        self.cylinders = []
        self.scales = []
        self.end = None
        return

    def checkDone(self):
        if self.end == None:
            raise ValueError('Geometry input after an END')

    def addWire(self, wobj):
        self.checkDone()
        if wobj.tag in self.wires:
            raise ValueError('Duplicate wire tag: <%s>' % wobj.line)
        self.wires[wobj.tag] = wobj
        self.wiresInOrder.append(wobj)
        self.lastWireAdded = wobj
        return

    def addHelix(self, hobj):
        self.checkDone()
        if hobj.tag in self.wires:
            raise ValueError('Duplicate wire tag: <%s>' % wobj.line)
        self.wires[wobj.tag] = wobj
        self.wiresInOrder.append(wobj)
        self.lastWireAdded = wobj
        return

    def addCylinder(self, cobj):
        self.checkDone()
        self.cylinders.append(cobj)
        # might need to process here and add additional wires, this one is complex
        return

    def addTransform(self, tobj):
        self.checkDone()
        self.transforms.append(wobj)
        # -spc- TODO I ~think that at this point we need to run each transform against the wires we
        # have to this point, in case tranforms and new wires are interleaved. Check xnec2c handling
        return

    def addScale(self, sobj):
        self.checkDone()
        self.scales.append(sobj)
        # -spc- TODO I ~think that at this point we need to run the scale against the wires we
        # have to this point, in case scales and new wires are interleaved. Check xnec2c handling
        return

    def addTaper(self, tobj):
        self.checkDone()
        self.lastWireAdded.addTaper(tobj)
        return

    def end(self, eobj):
        self.end = eobj
        return

class Excitation(object):
    """
    Excitation information for the model
    """
    def __init(self):
        return

# The overall model
class NecModel(object):
    """
    A class to encapsulate the information in an NEC2 antenna analysis model
    """

    def __init__(self, verbose=False):
        self.ops = {'CM':self.CM, 'CE':self.CE, 'GA':self.GA, 'GM':self.GM, 'GR':self.GR, 'GW':self.GW, 'GE':self.GE,
                    'GF':self.GF, 'GH':self.GH, 'GS':self.GS, 'GX':self.GX, 'SP':self.SP, 'SM':self.SM, 'SC':self.SC,
                    'GC':self.GC, 'CP':self.CP, 'EK':self.EK, 'EX':self.EX, 'FR':self.FR, 'GN':self.GN, 'KH':self.KH,
                    'LD':self.LD, 'GD':self.GD, 'NH':self.NH, 'NE':self.NE, 'NT':self.NT, 'NX':self.NX, 'RP':self.RP,
                    'PQ':self.PQ, 'PT':self.PT, 'TL':self.TL, 'WG':self.WG, 'XQ':self.XQ, 'EN':self.EN}
        if verbose:
            self.opts['verbose'] = True

        self.comments = []
        self.geometry = Geometry()

        return
    
    def inches(self, meters):
        return float(meters* 39.3700787)

    def CM(self, line):
        # Comment
        sline = line.strip()
        comm = sline[3:]
        if comm != '':
            self.comments.append(comm)
        if 'verbose' in self.opts:
            print 'Comment: %s' % comm
        return

    def CE(self, line):
        # End Comment
        sline = line.strip()
        comm = sline[3:]
        if comm != '':
            self.comments.append(comm)
        if verbose in self.opts:
            print 'Com End: %s' % comm
        return

    #
    # Structure geometry input
    #

    def GA(self, line):
        # GA Wire arc specification
        na = Arc(line)

        if 'verbose' in self.opts:
            na.dump()
        self.geometry.addWire(na)
        return

    def GM(self, line):
        # GM Coordinate transformation
        nt = Transform(line)
        if 'verbose' in self.opts:
            nt.dump()
        self.geometry.addTransform(nt)
        return

    def GR(self, line):
        # GR Generate cylindrical structure
        nc = Cylinder(line)
        if 'verbose' in self.opts:
            nc.dump()
        self.geometry.addCylinder(nc)
        return

    def GW(self, line):
        # GW Wire Specification
        nw = Wire(line)
        if 'verbose' in self.opts:
            nw.dump()
        self.geometry.addWire(nc)
        return

    def GC(self, line):
        # GC Wire Taper
        to = Taper(line)
        if 'verbose' in self.opts:
            to.dump()
        self.geometry.addTaper(to)
        return

    def GE(self, line):
        # GE End of geometry Input
        ge = GeometryEnd(line)
        if 'verbose' in self.opts:
            ge.dump()
        self.geometry.end(ge)
        return

    def GF(self, line):
        # GF Read Numerical Green's Function File
        raise ValueError("GF Read Numerical Green's Function File (not supported in xnec2c)")

    def GH(self, line):
        # GH Helix/Spiral specification
        nh = Helix(line)
        if 'verbose' in self.opts:
            nh.dump()
        self.geometry.addHelix(nh)
        return

    def GS(self, line):
        # GS Scale structure definitions
        ns = Scale(line)
        if 'verbose' in self.opts:
            ns.dump()
        self.geometry.addScale(ns)
        return

# zzzz End of geometry stuff

    def GX(self, line):
        # GX Reflection in coordinate plane
        parts = line.split()
        TNI = int(parts[1])
        REFL =  int(parts[2])
        # nec2c just assumes that the space-delimited number is extended
        # to three digits, and any non-zero digit is set to one
        dstr = '%03d' % REFL
        if 'verbose' in self.opts:
            print 'GX Reflection in coordinate plane:'
            print '    Tag Increment:      %d' % TNI
            if dstr[0] == '1':
                print '    Reflected along X axis'
            if dstr[1] == '1':
                print '    Reflected along Y axis'
            if dstr[2] == '1':
                print '    Reflected along Z axis'
        return

    def SP(self, line):
        # Surface Patch
        parts = line.split()
        NOP1 = int(parts[1])
        NS = int(parts[2])
        X1 = float(parts[3])
        Y1 = float(parts[4])
        Z1 = float(parts[5])
        X2 = float(parts[6])
        Y2 = float(parts[7])
        Z2 = float(parts[8])
        if 'verbose' in self.opts:
            print 'SP Surface Patch:'
            if NS == 0:
                # Arbitrary shape
                print '    Shape: Arbitrary'
            elif NS == 1:
                print '    Shape: Rectangular'
            elif NS == 2:
                print '    Shape: Triangular'
            elif NS == 3:
                print '    Shape: Quadrilateral'
            print '    X1     %f m (%.2f in)' % (X1, inches(X1))
            print '    Y1     %f m (%.2f in)' % (Y1, inches(Y1))
            print '    Z1     %f m (%.2f in)' % (Z1, inches(Z1))
            print '    X2     %f m (%.2f in)' % (X2, inches(X2))
            print '    Y2     %f m (%.2f in)' % (Y2, inches(Y2))
            print '    Z2     %f m (%.2f in)' % (Z2, inches(Z2))
        return

    def SM(self, line):
        # Multiple Patch Surface
        # Not tested
        parts = line.split()
        NX = int(parts[1])
        NY = int(parts[2])
        X1 = float(parts[3])
        Y1 = float(parts[4])
        Z1 = float(parts[5])
        X2 = float(parts[6])
        Y2 = float(parts[7])
        Z2 = float(parts[8])
        if 'verbose' in self.opts:
            print 'SM Multiple Patch Surface:'
            print '    Patches across X:     %d' % NX
            print '    Patches across Y:     %d' % NY
            print '    Corner 1 x:           %f m (%.2f in)' % (X1, inches(X1))
            print '    Corner 1 y:           %f m (%.2f in)' % (Y1, inches(Y1))
            print '    Corner 1 z:           %f m (%.2f in)' % (Z1, inches(Z1))
            print '    Corner 2 x:           %f m (%.2f in)' % (X2, inches(X2))
            print '    Corner 2 y:           %f m (%.2f in)' % (Y2, inches(Y2))
            print '    Corner 2 z:           %f m (%.2f in)' % (Z2, inches(Z2))
        return

    def SC(self, line):
        # Surface Multiple Patch Continuation
        # Not tested
        parts = line.split()
        NOP1 = int(parts[1])
        NS = int(parts[2])
        X3 = float(parts[3])
        Y3 = float(parts[4])
        Z3 = float(parts[5])
        if 'verbose' in self.opts:
            print 'SC Surface Multiple Patch Continuation:'
            print '    Corner 3 x:           %f m (%.2f in)' % (X3, inches(X3))
            print '    Corner 3 y:           %f m (%.2f in)' % (Y3, inches(Y3))
            print '    Corner 3 z:           %f m (%.2f in)' % (Z3, inches(Z3))
        return

    #
    # Program Control
    #

    def CP(self, line):
        # Maximum Coupling Calculation
        if 'verbose' in self.opts:
            print 'CP Maximum Coupling Calculation:'
            print '    Not supported in nec2c'
        return

    def EK(self, line):
        # Extended thin-wire kernel
        parts = line.split()
        I1 = int(parts[1])
        if 'verbose' in self.opts:
        print 'EK Extended thin-wire kernel:'
            if I1 == 0:
                print '    Initiate extended thin-wire kernel'
            elif I1 == -1:
                print '    Return to standard thin-wire kernel'
            else:
                print '    Unexpected value %d' % I1
        return

    def EX(self, line):
        # Excitation
        parts = line.split()
        TYPE = int(parts[1])
        I2   = int(parts[2])
        I3   = int(parts[3])
        I4   = int(parts[4])
        F1 = float(parts[5])
        F2 = float(parts[6])
        F3 = float(parts[7])
        if type not in [0,1,2,3,4,5]:
            if 'verbose' in self.opts:
                print '    Unknown Excitation type %d' % TYPE
            return
        if 'verbose' in self.opts:
            print 'EX Excitation:'
            if TYPE == 0:
                print '    Voltage Source'
            elif TYPE == 1:
                print '    Incident Plane Wave, linear polarization'
            elif TYPE == 2:
                print '    Incident Plane Wave, right hand elliptic polarization'
            elif TYPE == 3:
                print '    Incident Plane Wave, left hand elliptic polarization'
            elif TYPE == 4:
                print '    Elementary current source'
            elif TYPE == 5:
                print '    Voltage Source (current-slope-discontinuity)'

            if TYPE in [0,5]:
                # voltage
                print '    Driven Tag:      %d' % I2
                print '    Driven Segment:  %d' % I3
                # expand I4 to two positions
                tmp = '%02d' % I4
                if tmp[0] == '0':
                    col19 = False
                elif tmp[0] == '1':
                    col19 = True
                    print '    maximum relative admittance matrix asymmetry for source segment and network connection will be calculated and printed'
                else:
                    col19 = False
                    print 'Unexpected digit <%d> in field I4: %s' % (I4, line)
                if tmp[1] == '0':
                    col20 = False
                elif tmp[1] == '1':
                    col20 = True
                    print '    the input impedance at voltage sources is always printed directly before the segment currents in the output'
                else:
                    col20 = False
                    print 'Unexpected digit <%d> in field I4: %s' % (I4, line)
                print '    Real part:       %f volts' % F1
                print '    Imaginary part:  %f volts' % F2
                if (col20):
                    print '    Normalization constant for impedance:  %f' % F3

            elif TYPE in [1,2,3]:
                # plane wave
                print '    Number of Theta angles:      %d' % I2
                print '    Number of Phi angles:        %d' % I3
                tmp = '%01d' % I4
                if tmp[0] == '0':
                    pass
                elif tmp[0] == '1':
                    print '    maximum relative admittance matrix asymmetry for network connections will be calculated'
                else:
                    print 'Unexpected digit <%d> in field I4: %s' % (I4, line)
                print '    Incident Wave vector Theta:     %f degrees' % F1
                print '    Incident Wave vector Phi:       %f degrees' % F2
                print '    Polarization angle Eta:         %f degrees' % F3
                print '    Theta angle stepping increment: %f degrees' % F4
                print '    Phi angle stepping increment:   %f degrees' % F5
                print '    Ratio of minor to major Axes:   %f degrees' % F6
            else:
                # must be 4 - current source
                tmp = '%01d' % I4
                if tmp[0] == '0':
                    pass
                elif tmp[0] == '1':
                    print '    maximum relative admittance matrix asymmetry for network connections will be calculated'
                else:
                    print 'Unexpected digit <%d> in field I4: %s' % (I4, line)
                print '    X position:     %f m (%.2f in)' % (F1, inches(F1))
                print '    Y position:     %f m (%.2f in)' % (F2, inches(F2))
                print '    Z position:     %f m (%.2f in)' % (F3, inches(F3))
                print '    Alpha:          %f degrees'   % F4
                print '    Beta:           %f degrees'   % F5
                print '    Current Moment: %f' % F5
        return

    def FR(self, line):
        # Frequency
        parts = line.split()
        IFRQ = int(parts[1])
        NFRQ = int(parts[2])
        NOP1 = int(parts[3])
        NOP2 = int(parts[4])
        FREQ = float(parts[5])
        DELFRQ = float(parts[6])
        if 'verbose' in self.opts:
            print 'FR Frequency:'
            if IFRQ == 0:
                    print '    Linear Stepping'
                    endfreq = FREQ + (NFRQ * DELFRQ)
            elif IFRQ == 1:
                    print '    Multiplicative Stepping'
                    endfreq = FREQ
                    for i in range(NFRQ):
                        endfreq = endfreq * DELFRQ
            else:
                    print '    Unexpected stepping value %d: %s' % (IFRQ, line)
            print '    Number of Steps:    %d' % NFRQ
            print '    Frequency:          %f MHz'   % FREQ
            print '    Increment:          %f' % DELFRQ
            print '    Ending Frequency:   %f MHz [calculated]' % endfreq
        return

    def GN(self, line):
        # Ground parameters
        parts  = line.split()
        IPERF  = int(parts[1])
        if 'verbose' in self.opts:
            print 'GN Ground parameters:'
        if IPERF == -1:
            if 'verbose' in self.opts:
                print '    Free space (nullifies any previous ground definitions)'
            return
        NRADL  = int(parts[2])
        NOP1   = int(parts[3])
        NOP2   = int(parts[4])
        EPSE   = float(parts[5])
        SIG    = float(parts[6])
        if 'verbose' in self.opts:
            if IPERF == 0:
                    print '    Finite ground, reflection-coefficient approximation'
            elif IPERF == 1:
                    print '    Perfectly conducting ground'
            elif IPERF == 2:
                    print '    Finite ground, Sommerfield/Norton method'
            else:
                    print '    Unexpected ground-type flag %d: %s' % (IPERF, line)
            print '    Relative Dialectric Constant:  %f' % EPSE
            print '    Conductivity:                  %f (mhos/meter)' % SIG
            print '    Number of Radials:             %d' % NRADL
            if NRADL == 0:
                # cliff problem
                EPSR2   = float(parts[7])
                SIG2    = float(parts[8])
                CLT     = float(parts[9])
                CHT     = float(parts[10])
                print '    Relative Dialectric Constant (second medium):   %f'   % EPSR2
                print '    Conductivity (second medium):                   %f mhos/meter'   % SIG2
                print '    Distance from origin to join between mediums:   %f m (%.2f in)' % (CLT, inches(CLT))
                print '    Distance medium 2 is below medium 1:            %f m (%.2f in)' % (CHT, inches(CHT))
            else:
                # we have radials
                SRAD   = float(parts[7])
                WRAD    = float(parts[8])
                print '    Radius of radial screen:       %f m (%.2f in)' % (SRAD, inches(SRAD))
                print '    Radius of radial wires:        %f m (%.2f in)' % (WRAD, inches(WRAD))
        return
zzzzzzz
    def KH(self, line):
        # Interaction Approximation Range
        print 'KH Interaction Approximation Range:'
        #print '    %s' % line,
        parts  = line.split()
        RKH  = int(parts[1])
        print '    Approximation used for interactions over %d wavelengths' % RKH
        return

    def LD(self, line):
        # Loading
        print 'LD Loading:'
        #print '    %s' % line,
        parts = line.split()
        LDTYP = int(parts[1])
        if LDTYP == -1:
                print '    Short all loads (nullifies any previous loads)'
                return
        if LDTYP == 0:
                print '    Series RLC, input Ohms, Henries, Farads'
        elif LDTYP == 1:
                print '    Parallel RLC, input Ohms, Henries, Farads'
        elif LDTYP == 2:
                print '    Series RLC, input Ohms/meter, Henries/meter, Farads/meter'
        elif LDTYP == 3:
                print '    Parallel RLC, input Ohms/meter, Henries/meter, Farads/meter'
        elif LDTYP == 4:
                print '    Impedance, input resistance and reactance in Ohms'
        elif LDTYP == 5:
                print '    Wire conductivity, mhos/meter'
        else:
                print '    Unexpected value of %d for Load Type: %s' % (LDTYP, line)

        LDTAG = int(parts[2])
        LDTAGF = int(parts[3])
        LDTAGT = int(parts[4])
        print '    Tag:             %d' % LDTAG
        print '    Start segment:   %d' % LDTAGF
        print '    End segment:     %d' % LDTAGT

        ZLR = float(parts[5])
        ZLI = float(parts[6])
        ZLC = float(parts[7])
        if LDTYP in [0,1]:
            units1 = 'Ohms'
            units2 = 'Henries'
            units3 = 'Farads'
            print '    Resistance:  %f %s' % (ZLR, units1)
            print '    Inductance:  %f %s' % (ZLI, units2)
            print '    Capacitance: %f %s' % (ZLC, units3)
        elif LDTYP in [2,3]:
            units1 = 'Ohms/meter'
            units2 = 'Henries/meter'
            units3 = 'Farads/meter'
            print '    Resistance:  %f %s' % (ZLR, units1)
            print '    Inductance:  %f %s' % (ZLI, units2)
            print '    Capacitance: %f %s' % (ZLC, units3)
        elif LDTYP == 4:
            ohms = 'Ohms'
            print '    Resistance:  %f %s' % (ZLR, ohms)
            print '    Reactance:   %f %s' % (ZLI, ohms)
        else:
            # LDTYP=5
            mhosmeter = 'mhos/meter'
            print '    Conductivity:  %f %s' % (ZLR, mhosmeter)
        return

    def GD(self, line):
        # Additional ground parameters
        print 'GD Additional ground parameters:'
        #print '    %s' % line,
        parts  = line.split()
        NOP1   = int(parts[1])
        NOP2   = int(parts[2])
        NOP3   = int(parts[3])
        NOP4   = int(parts[4])
        EPSR2  = float(parts[5])
        SIG2   = float(parts[6])
        CLT    = float(parts[7])
        CHT    = float(parts[8])
        print '    Relative Dialectric Constant (second medium):   %f'   % EPSR2
        print '    Conductivity (second medium):                   %f mhos/meter'   % SIG2
        print '    Distance from origin to join between mediums:   %f m (%.2f in)' % (CLT, inches(CLT))
        print '    Distance medium 2 is below medium 1:            %f m (%.2f in)' % (CHT, inches(CHT))
        return

    def nhne(self, line):
        #print '    %s' % line,
        parts  = line.split()
        NEAR   = int(parts[1])
        NRX    = int(parts[2])
        NRY    = int(parts[3])
        NRZ    = int(parts[4])
        XNR    = float(parts[5])
        YNR    = float(parts[6])
        ZNR    = float(parts[7])
        DXNR   = float(parts[8])
        DYNR   = float(parts[9])
        DZNR   = float(parts[10])

        if NEAR == 0:
            print '    Using Rectangular Coordinates'
            print '    Points desired in X:     %d'   % NRX
            print '    Points desired in Y:     %d'   % NRY
            print '    Points desired in Z:     %d'   % NRZ
            print '    First field point X:     %f m (%.2f in)' % (XNR, inches(XNR))
            print '    First field point Y:     %f m (%.2f in)' % (YNR, inches(YNR))
            print '    First field point Z:     %f m (%.2f in)' % (ZNR, inches(ZNR))
            print '    Stepping increment in X  %f m (%.2f in)' % (DXNR, inches(DXNR))
            print '    Stepping increment in Y  %f m (%.2f in)' % (DYNR, inches(DYNR))
            print '    Stepping increment in Z  %f m (%.2f in)' % (DZNR, inches(DZNR))

        elif NEAR == 1:
            print '    Using Spherical Coordinates'
            print '    Points desired in r:         %d'   % NRX
            print '    Points desired in phi:       %d'   % NRY
            print '    Points desired in theta:     %d'   % NRZ
            print '    First field point r:         %f m (%.2f in)' % (XNR, inches(XNR))
            print '    First field point phi:       %f degrees' % YNR
            print '    First field point theta:     %f degrees' % ZNR
            print '    Stepping increment in r:     %f m (%.2f in)' % (DXNR, inches(DXNR))
            print '    Stepping increment in phi:   %f degrees' % DYNR
            print '    Stepping increment in theta: %f degrees' % DZNR
        else:
            print '    Unexpected coordinate type argument %d: %s' % (NEAR, line)
        return

    def NH(self, line):
        # Near Field NH
        print 'NH Near Magneticic Field:'
        nhne(self, line)
        return

    def NE(self, line):
        # Near Field NE
        print 'NE Near Electric Field:'
        nhne(self, line)
        return

    def NT(self, line):
        # Networks
        # Untested
        print 'NT Networks:'
        #print '    %s' % line,
        parts  = line.split()
        TAG1   = int(parts[1])
        SEG1   = int(parts[2])
        TAG2   = int(parts[3])
        SEG2   = int(parts[4])
        Y11R   = float(parts[5])
        Y11I   = float(parts[6])
        Y12R   = float(parts[7])
        Y12I   = float(parts[8])
        Y22R   = float(parts[9])
        Y22I   = float(parts[10])
        print '    Port 1 tag:       %d' % TAG1
        print '    Port 1 segment:   %d' % SEG1
        print '    Port 2 tag:       %d' % TAG2
        print '    Port 2 segment:   %d' % SEG2
        print '    Real part of element (1, 1):        %f mhos' % Y11R
        print '    Imaginary part of element (1, 1):   %f mhos' % Y11I
        print '    Real part of element (1, 2):        %f mhos' % Y12R
        print '    Imaginary part of element (1, 2):   %f mhos' % Y12I
        print '    Real part of element (2, 2):        %f mhos' % Y22R
        print '    Imaginary part of element (2, 2):   %f mhos' % Y22I
        return

    def NX(self, line):
        # Next Structure
        print 'NX Next Structure:'
        #print '    %s' % line,
        print '    Not supported in nec2c'
        return

    def RP(self, line):
        # Radiation Pattern
        print 'RP Radiation Pattern:'
        #print '    %s' % line,
        parts  = line.split()
        MODE   = int(parts[1])
        if MODE == 0:
            print '    Mode: Normal (space-wave fields)'
        elif MODE == 1:
            print '    Mode: Surface wave propagating along ground'
        elif MODE == 2:
            print '    Mode: Linear cliff with antenna above upper level'
        elif MODE == 3:
            print '    Mode: Circular cliff centered at origin, antenna above upper level'
        elif MODE == 4:
            print '    Mode: Radial wire ground screen centered at origin'
        elif MODE == 5:
            print '    Mode: Radial wire ground screen and linear cliff'
        elif MODE == 6:
            print '    Mode: Radial wire ground screen and circular cliff'
        NTH    = int(parts[2])
        NPH    = int(parts[3])
        XNDA   = int(parts[4])
        THETS   = float(parts[5])
        PHIS   = float(parts[6])
        DTH   = float(parts[7])
        DPH   = float(parts[8])
        RFLD   = float(parts[9])
        GNOR   = float(parts[10])
        if MODE == 1:
            print '    Points in z:           %d'   % NTH
        else:
            print '    Points in theta:       %d'   % NTH
        print '    Points in phi:         %d'   % NPH
        xnda = '%04d' % XNDA
        if xnda[0] == '0':
            print '    major-axis, minor-axis, and total gain printed'
        elif xnda[0] == '1':
            print '    vertical, horizontal, and total gain printed'
        else:
            print '    Unexpected value in XNDA field %s : %s' % (XNDA, line)

        if xnda[1] == '0':
            print '    No normalized gain'
        elif xnda[1] == '1':
            print '    Major axis gain normalized'
        elif xnda[1] == '2':
            print '    Minor axis gain normalized'
        elif xnda[1] == '3':
            print '    Vertical axis gain normalized'
        elif xnda[1] == '4':
            print '    Horizontal axis gain normalized'
        elif xnda[1] == '5':
            print '    Total gain normalized'
        else:
            print '    Unexpected value in XNDA field %s : %s' % (XNDA, line)

        if xnda[2] == '0':
            print '    Use power Gain'
        elif xnda[2] == '1':
            print '    Use directive gain'
        else:
            print '    Unexpected value in XNDA field %s : %s' % (XNDA, line)

        if xnda[3] == '0':
            print '    No averaging'
        elif xnda[3] == '1':
            print '    Average gain computed'
        elif xnda[3] == '2':
            print '    Average gain computed, supress printing'
        else:
            print '    Unexpected value in XNDA field %s : %s' % (XNDA, line)

        if MODE == 1:
            print '    Initial Z coordinate:      %f m (%.2f in)' % (THETS, inches(THETS))
            print '    Initial phi angle:         %f degrees' % PHIS
            print '    Increment for z:           %f m (%.2f in)' % (DTH, inches(DTH))
            print '    Increment for phi:         %f degrees' % DPH
            print '    Cylindrical coord phi:     %f m (%.2f in)' % (RFLD, inches(RFLD))
            print '    Range of z:                %f - %f m, %.2f - %.2f (in) [calculated]' % (THETS, THETS + (NTH * DTH), inches(THETS), inches(THETS + (NTH * DTH)))
            print '    Range of phi:              %f - %f degrees [calculated]' % (PHIS, PHIS + (NPH * DPH))
        else:
            print '    Initial theta angle:       %f degrees' % THETS
            print '    Initial phi angle:         %f degrees' % PHIS
            print '    Increment for theta:       %f degrees' % DTH
            print '    Increment for phi:         %f degrees' % DPH
            print '    Radial distance:           %f m (%.2f in)' % (RFLD, inches(RFLD))
            print '    Range of theta:            %f - %f degrees [calculated]' % (THETS, THETS + (NTH * DTH))
            print '    Range of phi:              %f - %f degrees [calculated]' % (PHIS, PHIS + (NPH * DPH))
        print '    Gain normalization factor: %f' % GNOR
        return

    def PQ(self, line):
        # Print control for charge on wires
        print 'PQ Print control for charge on wires:'
        #print '    %s' % line,
        print '    Not supported in nec2c'
        return

    def PT(self, line):
        # Page Title / Print Control for Current on Wires 
        print 'PT Page Title / Print Control for Current on Wires :'
        #print '    %s' % line,
        print '    Not supported in nec2c'
        return

    def TL(self, line):
        # Transmission Line
        print 'TL Transmission Line:'
        #print '    %s' % line,
        parts  = line.split()
        TAG1   = int(parts[1])
        SEG1   = int(parts[2])
        TAG2   = int(parts[3])
        SEG2   = int(parts[4])

        IMP    = float(parts[5])
        LEN    = float(parts[6])
        R1     = float(parts[7])
        I1     = float(parts[8])
        R2     = float(parts[9])
        I2     = float(parts[10])
        print '    Port 1 tag:       %d' % TAG1
        print '    Port 1 segment:   %d' % SEG1
        print '    Port 2 tag:       %d' % TAG2
        print '    Port 2 segment:   %d' % SEG2

        print '    Characteristic impedance:                   %f ohms' % IMP
        print '    Length:                                     %f m (%.2f in)' % (LEN, inches(LEN))
        print '    Real part of shunt admittance, end 1        %f mhos' % R1
        print '    Imaginary part of shunt admittance, end 1   %f mhos' % I1
        print '    Real part of shunt admittance, end 2        %f mhos' % R2
        print '    Imaginary part of shunt admittance, end 2   %f mhos' % I2
        return

    def WG(self, line):
        # Write NGF file
        print 'WG Write NGF file:'
        #print '    %s' % line,
        print '    Not supported in nec2c'
        return

    def XQ(self, line):
        # Execute
        print 'XQ Execute:'
        #print '    %s' % line,
        print '    Ignored in nec2c'
        return

    def EN(self, line):
        # End of run
        print 'EN End of run'
        #print '    %s' % line,
        #print 'End of Run'
        return

def main():
    printops = []
    if len(sys.argv) == 1:
        print 'Usage: %s filename [taglist]' % sys.argv[0]
        print '   where filename is mandatory and [taglist] is an optional comma-delimited list of tags to print (no spaces in taglist)'
        return

    if len(sys.argv) > 2:
        opslist = sys.argv[2]
        for op in opslist.split(','):
            printops.append(op.upper())
    inf = open(sys.argv[1])
    for line in inf.readlines():
        lparts = line.split()
        op = lparts[0]
        if (len(printops) == 0) or (len(printops) and op in printops):
            if op  in ops:
                ops[lparts[0].upper()](line)
            else:
                print 'Unhandled -> ', line,
    return

main()


# vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
