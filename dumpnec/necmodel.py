# Comment #!/usr/bin/python3
import sys

def inches(meters):
    return float(meters * 39.3700787)

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

    def dump(self):
        print("GC Wire Taper")
        print('    Ratio to prior segment:   %f' % self.ratio)
        print('    First Segment Radius:      %f m (%.2f in)' % (self.radiusFirst, inches(self.radiusFirst)))
        print('    Last  Segment Radius:      %f m (%.2f in)' % (self.radiusLast, inches(self.radiusLast)))

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
        print('GW Wire Specification:')
        print('    Tag:           %d' % self.tag)
        print('    Segments:      %d' % self.segments)
        print('    Start x:       %f m (%.2f in)' % (self.startX, inches(self.startX)))
        print('    Start y:       %f m (%.2f in)' % (self.startY, inches(self.startY)))
        print('    Start z:       %f m (%.2f in)' % (self.startZ, inches(self.startZ)))
        print('    End   x:       %f m (%.2f in)' % (self.endX, inches(self.endX)))
        print('    End   y:       %f m (%.2f in)' % (self.endY, inches(self.endY)))
        print('    End   z:       %f m (%.2f in)' % (self.endZ, inches(self.endZ)))
        print('    Wire Radius:   %f m (%.2f in)' % (self.wireRadius, inches(self.wireRadius)))
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
        print('GA Wire arc specification:')
        print('    Tag:          %d' % self.tag)
        print('    Segments:     %d' % self.segments)
        print('    Arc Radius:   %f m (%.2f in)' % (self.radius, inches(self.radius)))
        print('    Angle start:  %f degrees' % self.angleStart)
        print('    Angle end:    %f degrees' % self.AngleEnd)
        print('    Wire Radius:  %f' % self.wireRadius)
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

    def dump(self):
        print('GH Helix/Spiral specification:')
        if self.length == 0:
            print('    (Spiral)')
        elif self.length > 0:
            print('    (Right Handed Helix)')
        else:
            print('    (Left Handed Helix)')
        print('    Segment Tag:         %d' % self.tag)
        print('    Segments:            %d' % self.segments)
        print('    Turn Spacing:        %f m (%.2f in)' % (self.spacing, inches(self.spacing)))
        print('    Helix Length HL:     %f m (%.2f in)' % (self.length, inches(self.length)))
        print('    Radius in x @ z=0:   %f m (%.2f in)' % (self.radiusXZ0, inches(self.radiusXZ0)))
        print('    Radius in y @ z=0:   %f m (%.2f in)' % (self.radiusYZ0, inches(self.radiusYZ0)))
        print('    Radius in x @ z=HL:  %f m (%.2f in)' % (self.radiusXZhl, inches(self.radiusXZhl)))
        print('    Radius in y @ z=HL:  %f m (%.2f in)' % (self.radiusYZhl, inches(self.radiusYZhl)))
        print('    Wire Radius:         %f' % self.radius)
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
            print('GR Generate cylindrical structure:')
            print('    Tag Increment:                     %d' % self.tagIncrement)
            print('    Number of times structure occurs:  %d' % self.numOccurrances)
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
        print('GM Coordinate transformation:')
        print('    Tag Increment:                %d' % self.tagIncrement)
        print('    New Structs                   %d' % self.newStructs)
        print('    Rotation about x:             %f degrees' % self.rotationX)
        print('    Rotation about y:             %f degrees' % self.rotationY)
        print('    Rotation about z:             %f degrees' % self.rotationZ)
        print('    Translation in x:             %f m (%.2f in)' % (self.translationX, inches(self.translationX)))
        print('    Translation in y:             %f m (%.2f in)' % (self.translationY, inches(self.translationY)))
        print('    Translation in z:             %f m (%.2f in)' % (self.translationZ, inches(self.translationZ)))
        print('    Initial Translation Segment:  %d' % int(round(self.its)))
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
        print('GS Scale Structure Domensions:')
        print('    Scale Factor:                 %f' % self.factor)
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

    def dump(self):
        print('    Ground Plane Flag:   %d' % self.groundPlane, end='')
        if self.groundPlane == 0:
            print(' (No Ground Plane Present)')
        elif self.groundPlane == 1:
            print(' (Ground Plane - Current Expansion Modified)')
        elif self.groundPlane == -1:
            print(' (Ground Plane - Current Expansion Not Modified)')
        else:
            print(' (Unknown Value)')
        return

# Classes for each subsection
class Geometry(object):
    """
    Antenna geometry information
    """
    def __init__(self):
        self.lastWireAdded = None
        self.wiresInOrder = []
        self.wires = {}
        self.transforms = []
        self.cylinders = []
        self.scales = []
        self.end = None
        return

    def checkDone(self):
        if self.end is not None:
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
        self.wires[hobj.tag] = hobj
        self.wiresInOrder.append(hobj)
        self.lastWireAdded = hobj
        return

    def addCylinder(self, cobj):
        self.checkDone()
        self.cylinders.append(cobj)
        # might need to process here and add additional wires, this one is complex
        return

    def addTransform(self, tobj):
        self.checkDone()
        self.transforms.append(tobj)
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

    def dump(self):
        for item in self.wiresInOrder:
            item.dump()
        for t in self.transforms:
            t.dump()
        for c in self.cylinders:
            c.dump()
        for s in self.scales:
            s.dump()
        if self.end is not None:
            self.end.dump()
        return

class Frequency(object):
    """
    Frequency sweep definition (FR card)
    """
    def __init__(self, line):
        parts = line.split()
        self.stepping  = int(parts[1])    # 0=linear, 1=multiplicative
        self.numSteps  = int(parts[2])
        self.frequency = float(parts[5])  # MHz
        self.increment = float(parts[6])
        return

    def dump(self):
        print('FR Frequency:')
        if self.stepping == 0:
            print('    Linear Stepping')
            endfreq = self.frequency + (self.numSteps * self.increment)
        elif self.stepping == 1:
            print('    Multiplicative Stepping')
            endfreq = self.frequency
            for _ in range(self.numSteps):
                endfreq = endfreq * self.increment
        else:
            print('    Unexpected stepping value %d' % self.stepping)
            endfreq = self.frequency
        print('    Number of Steps:    %d' % self.numSteps)
        print('    Frequency:          %f MHz' % self.frequency)
        print('    Increment:          %f' % self.increment)
        print('    Ending Frequency:   %f MHz [calculated]' % endfreq)
        return

class Excitation(object):
    """
    Excitation source (EX card)
    """
    def __init__(self, line):
        parts = line.split()
        self.type    = int(parts[1])
        self.tag     = int(parts[2])
        self.segment = int(parts[3])
        self.i4      = int(parts[4])
        self.f1 = float(parts[5])
        self.f2 = float(parts[6])
        self.f3 = float(parts[7])
        # plane wave and current source types use additional fields
        if self.type in [1, 2, 3, 4]:
            self.f4 = float(parts[8])  if len(parts) > 8  else 0.0
            self.f5 = float(parts[9])  if len(parts) > 9  else 0.0
        if self.type in [1, 2, 3]:
            self.f6 = float(parts[10]) if len(parts) > 10 else 0.0
        return

    def dump(self):
        print('EX Excitation:')
        if self.type == 0:
            print('    Voltage Source')
        elif self.type == 1:
            print('    Incident Plane Wave, linear polarization')
        elif self.type == 2:
            print('    Incident Plane Wave, right hand elliptic polarization')
        elif self.type == 3:
            print('    Incident Plane Wave, left hand elliptic polarization')
        elif self.type == 4:
            print('    Elementary current source')
        elif self.type == 5:
            print('    Voltage Source (current-slope-discontinuity)')
        else:
            print('    Unknown Excitation type %d' % self.type)
            return
        if self.type in [0, 5]:
            print('    Driven Tag:      %d' % self.tag)
            print('    Driven Segment:  %d' % self.segment)
            tmp = '%02d' % self.i4
            col19 = (tmp[0] == '1')
            col20 = (tmp[1] == '1')
            if col19:
                print('    maximum relative admittance matrix asymmetry for source segment and network connection will be calculated and printed')
            if col20:
                print('    the input impedance at voltage sources is always printed directly before the segment currents in the output')
            print('    Real part:       %f volts' % self.f1)
            print('    Imaginary part:  %f volts' % self.f2)
            if col20:
                print('    Normalization constant for impedance:  %f' % self.f3)
        elif self.type in [1, 2, 3]:
            tmp = '%01d' % self.i4
            if tmp[0] == '1':
                print('    maximum relative admittance matrix asymmetry for network connections will be calculated')
            print('    Number of Theta angles:         %d' % self.tag)
            print('    Number of Phi angles:           %d' % self.segment)
            print('    Incident Wave vector Theta:     %f degrees' % self.f1)
            print('    Incident Wave vector Phi:       %f degrees' % self.f2)
            print('    Polarization angle Eta:         %f degrees' % self.f3)
            print('    Theta angle stepping increment: %f degrees' % self.f4)
            print('    Phi angle stepping increment:   %f degrees' % self.f5)
            print('    Ratio of minor to major Axes:   %f degrees' % self.f6)
        else:
            # type == 4
            tmp = '%01d' % self.i4
            if tmp[0] == '1':
                print('    maximum relative admittance matrix asymmetry for network connections will be calculated')
            print('    X position:     %f m (%.2f in)' % (self.f1, inches(self.f1)))
            print('    Y position:     %f m (%.2f in)' % (self.f2, inches(self.f2)))
            print('    Z position:     %f m (%.2f in)' % (self.f3, inches(self.f3)))
            print('    Alpha:          %f degrees' % self.f4)
            print('    Beta:           %f degrees' % self.f5)
            print('    Current Moment: %f' % self.f5)
        return

class Ground(object):
    """
    Ground parameters (GN card)
    """
    def __init__(self, line):
        parts = line.split()
        self.type = int(parts[1])  # -1=free space, 0=reflect coeff, 1=perfect, 2=Sommerfeld
        self.numRadials = 0
        self.dielectric = 0.0
        self.conductivity = 0.0
        self.radialRadius = 0.0
        self.radialWireRadius = 0.0
        self.cliffDielectric = 0.0
        self.cliffConductivity = 0.0
        self.cliffDist = 0.0
        self.cliffHeight = 0.0
        if self.type == -1:
            return
        self.numRadials   = int(parts[2])
        self.dielectric   = float(parts[5])
        self.conductivity = float(parts[6])
        if self.numRadials == 0:
            # cliff problem
            self.cliffDielectric   = float(parts[7])
            self.cliffConductivity = float(parts[8])
            self.cliffDist         = float(parts[9])
            self.cliffHeight       = float(parts[10])
        else:
            self.radialRadius     = float(parts[7])
            self.radialWireRadius = float(parts[8])
        return

    def dump(self):
        print('GN Ground parameters:')
        if self.type == -1:
            print('    Free space (nullifies any previous ground definitions)')
            return
        if self.type == 0:
            print('    Finite ground, reflection-coefficient approximation')
        elif self.type == 1:
            print('    Perfectly conducting ground')
        elif self.type == 2:
            print('    Finite ground, Sommerfeld/Norton method')
        else:
            print('    Unexpected ground-type flag %d' % self.type)
        print('    Relative Dielectric Constant:  %f' % self.dielectric)
        print('    Conductivity:                  %f (mhos/meter)' % self.conductivity)
        print('    Number of Radials:             %d' % self.numRadials)
        if self.numRadials == 0:
            print('    Relative Dielectric Constant (second medium):   %f' % self.cliffDielectric)
            print('    Conductivity (second medium):                   %f mhos/meter' % self.cliffConductivity)
            print('    Distance from origin to join between mediums:   %f m (%.2f in)' % (self.cliffDist, inches(self.cliffDist)))
            print('    Distance medium 2 is below medium 1:            %f m (%.2f in)' % (self.cliffHeight, inches(self.cliffHeight)))
        else:
            print('    Radius of radial screen:       %f m (%.2f in)' % (self.radialRadius, inches(self.radialRadius)))
            print('    Radius of radial wires:        %f m (%.2f in)' % (self.radialWireRadius, inches(self.radialWireRadius)))
        return

class AdditionalGround(object):
    """
    Additional ground parameters (GD card)
    """
    def __init__(self, line):
        parts = line.split()
        self.dielectric   = float(parts[5])
        self.conductivity = float(parts[6])
        self.cliffDist    = float(parts[7])
        self.cliffHeight  = float(parts[8])
        return

    def dump(self):
        print('GD Additional ground parameters:')
        print('    Relative Dielectric Constant (second medium):   %f' % self.dielectric)
        print('    Conductivity (second medium):                   %f mhos/meter' % self.conductivity)
        print('    Distance from origin to join between mediums:   %f m (%.2f in)' % (self.cliffDist, inches(self.cliffDist)))
        print('    Distance medium 2 is below medium 1:            %f m (%.2f in)' % (self.cliffHeight, inches(self.cliffHeight)))
        return

class Loading(object):
    """
    Segment loading (LD card)
    """
    def __init__(self, line):
        parts = line.split()
        self.type = int(parts[1])  # -1=clear all, 0-5=type
        if self.type == -1:
            return
        self.tag      = int(parts[2])
        self.segFirst = int(parts[3])
        self.segLast  = int(parts[4])
        self.z1 = float(parts[5])
        self.z2 = float(parts[6])
        self.z3 = float(parts[7])
        return

    def dump(self):
        print('LD Loading:')
        if self.type == -1:
            print('    Short all loads (nullifies any previous loads)')
            return
        if self.type == 0:
            print('    Series RLC, input Ohms, Henries, Farads')
        elif self.type == 1:
            print('    Parallel RLC, input Ohms, Henries, Farads')
        elif self.type == 2:
            print('    Series RLC, input Ohms/meter, Henries/meter, Farads/meter')
        elif self.type == 3:
            print('    Parallel RLC, input Ohms/meter, Henries/meter, Farads/meter')
        elif self.type == 4:
            print('    Impedance, input resistance and reactance in Ohms')
        elif self.type == 5:
            print('    Wire conductivity, mhos/meter')
        else:
            print('    Unexpected value of %d for Load Type' % self.type)
        print('    Tag:             %d' % self.tag)
        print('    Start segment:   %d' % self.segFirst)
        print('    End segment:     %d' % self.segLast)
        if self.type in [0, 1]:
            print('    Resistance:  %f Ohms' % self.z1)
            print('    Inductance:  %f Henries' % self.z2)
            print('    Capacitance: %f Farads' % self.z3)
        elif self.type in [2, 3]:
            print('    Resistance:  %f Ohms/meter' % self.z1)
            print('    Inductance:  %f Henries/meter' % self.z2)
            print('    Capacitance: %f Farads/meter' % self.z3)
        elif self.type == 4:
            print('    Resistance:  %f Ohms' % self.z1)
            print('    Reactance:   %f Ohms' % self.z2)
        else:
            print('    Conductivity:  %f mhos/meter' % self.z1)
        return

class Network(object):
    """
    Two-port network (NT card)
    """
    def __init__(self, line):
        parts = line.split()
        self.tag1 = int(parts[1])
        self.seg1 = int(parts[2])
        self.tag2 = int(parts[3])
        self.seg2 = int(parts[4])
        self.y11r = float(parts[5])
        self.y11i = float(parts[6])
        self.y12r = float(parts[7])
        self.y12i = float(parts[8])
        self.y22r = float(parts[9])
        self.y22i = float(parts[10])
        return

    def dump(self):
        print('NT Networks:')
        print('    Port 1 tag:       %d' % self.tag1)
        print('    Port 1 segment:   %d' % self.seg1)
        print('    Port 2 tag:       %d' % self.tag2)
        print('    Port 2 segment:   %d' % self.seg2)
        print('    Real part of element (1, 1):        %f mhos' % self.y11r)
        print('    Imaginary part of element (1, 1):   %f mhos' % self.y11i)
        print('    Real part of element (1, 2):        %f mhos' % self.y12r)
        print('    Imaginary part of element (1, 2):   %f mhos' % self.y12i)
        print('    Real part of element (2, 2):        %f mhos' % self.y22r)
        print('    Imaginary part of element (2, 2):   %f mhos' % self.y22i)
        return

class TransmissionLine(object):
    """
    Transmission line (TL card)
    """
    def __init__(self, line):
        parts = line.split()
        self.tag1 = int(parts[1])
        self.seg1 = int(parts[2])
        self.tag2 = int(parts[3])
        self.seg2 = int(parts[4])
        self.impedance = float(parts[5])
        self.length    = float(parts[6])
        self.shunt1r   = float(parts[7])
        self.shunt1i   = float(parts[8])
        self.shunt2r   = float(parts[9])
        self.shunt2i   = float(parts[10])
        return

    def dump(self):
        print('TL Transmission Line:')
        print('    Port 1 tag:       %d' % self.tag1)
        print('    Port 1 segment:   %d' % self.seg1)
        print('    Port 2 tag:       %d' % self.tag2)
        print('    Port 2 segment:   %d' % self.seg2)
        print('    Characteristic impedance:                   %f ohms' % self.impedance)
        print('    Length:                                     %f m (%.2f in)' % (self.length, inches(self.length)))
        print('    Real part of shunt admittance, end 1        %f mhos' % self.shunt1r)
        print('    Imaginary part of shunt admittance, end 1   %f mhos' % self.shunt1i)
        print('    Real part of shunt admittance, end 2        %f mhos' % self.shunt2r)
        print('    Imaginary part of shunt admittance, end 2   %f mhos' % self.shunt2i)
        return

class RadiationPattern(object):
    """
    Radiation pattern request (RP card)
    """
    def __init__(self, line):
        parts = line.split()
        self.mode       = int(parts[1])
        self.numTheta   = int(parts[2])
        self.numPhi     = int(parts[3])
        self.xnda       = int(parts[4])
        self.thetaStart = float(parts[5])
        self.phiStart   = float(parts[6])
        self.thetaStep  = float(parts[7])
        self.phiStep    = float(parts[8])
        self.radialDist = float(parts[9])
        self.gainNorm   = float(parts[10])
        return

    def dump(self):
        print('RP Radiation Pattern:')
        if self.mode == 0:
            print('    Mode: Normal (space-wave fields)')
        elif self.mode == 1:
            print('    Mode: Surface wave propagating along ground')
        elif self.mode == 2:
            print('    Mode: Linear cliff with antenna above upper level')
        elif self.mode == 3:
            print('    Mode: Circular cliff centered at origin, antenna above upper level')
        elif self.mode == 4:
            print('    Mode: Radial wire ground screen centered at origin')
        elif self.mode == 5:
            print('    Mode: Radial wire ground screen and linear cliff')
        elif self.mode == 6:
            print('    Mode: Radial wire ground screen and circular cliff')
        if self.mode == 1:
            print('    Points in z:           %d' % self.numTheta)
        else:
            print('    Points in theta:       %d' % self.numTheta)
        print('    Points in phi:         %d' % self.numPhi)
        xnda = '%04d' % self.xnda
        if xnda[0] == '0':
            print('    major-axis, minor-axis, and total gain printed')
        elif xnda[0] == '1':
            print('    vertical, horizontal, and total gain printed')
        else:
            print('    Unexpected value in XNDA field %s' % self.xnda)
        if xnda[1] == '0':
            print('    No normalized gain')
        elif xnda[1] == '1':
            print('    Major axis gain normalized')
        elif xnda[1] == '2':
            print('    Minor axis gain normalized')
        elif xnda[1] == '3':
            print('    Vertical axis gain normalized')
        elif xnda[1] == '4':
            print('    Horizontal axis gain normalized')
        elif xnda[1] == '5':
            print('    Total gain normalized')
        else:
            print('    Unexpected value in XNDA field %s' % self.xnda)
        if xnda[2] == '0':
            print('    Use power gain')
        elif xnda[2] == '1':
            print('    Use directive gain')
        else:
            print('    Unexpected value in XNDA field %s' % self.xnda)
        if xnda[3] == '0':
            print('    No averaging')
        elif xnda[3] == '1':
            print('    Average gain computed')
        elif xnda[3] == '2':
            print('    Average gain computed, suppress printing')
        else:
            print('    Unexpected value in XNDA field %s' % self.xnda)
        if self.mode == 1:
            print('    Initial Z coordinate:      %f m (%.2f in)' % (self.thetaStart, inches(self.thetaStart)))
            print('    Initial phi angle:         %f degrees' % self.phiStart)
            print('    Increment for z:           %f m (%.2f in)' % (self.thetaStep, inches(self.thetaStep)))
            print('    Increment for phi:         %f degrees' % self.phiStep)
            print('    Cylindrical coord phi:     %f m (%.2f in)' % (self.radialDist, inches(self.radialDist)))
            print('    Range of z:                %f - %f m, %.2f - %.2f in [calculated]' % (self.thetaStart, self.thetaStart + (self.numTheta * self.thetaStep), inches(self.thetaStart), inches(self.thetaStart + (self.numTheta * self.thetaStep))))
            print('    Range of phi:              %f - %f degrees [calculated]' % (self.phiStart, self.phiStart + (self.numPhi * self.phiStep)))
        else:
            print('    Initial theta angle:       %f degrees' % self.thetaStart)
            print('    Initial phi angle:         %f degrees' % self.phiStart)
            print('    Increment for theta:       %f degrees' % self.thetaStep)
            print('    Increment for phi:         %f degrees' % self.phiStep)
            print('    Radial distance:           %f m (%.2f in)' % (self.radialDist, inches(self.radialDist)))
            print('    Range of theta:            %f - %f degrees [calculated]' % (self.thetaStart, self.thetaStart + (self.numTheta * self.thetaStep)))
            print('    Range of phi:              %f - %f degrees [calculated]' % (self.phiStart, self.phiStart + (self.numPhi * self.phiStep)))
        print('    Gain normalization factor: %f' % self.gainNorm)
        return

class NearField(object):
    """
    Near field request (NE or NH card)
    """
    def __init__(self, line, fieldType):
        parts = line.split()
        self.fieldType  = fieldType     # 'E' or 'H'
        self.coordType  = int(parts[1]) # 0=rectangular, 1=spherical
        self.nx = int(parts[2])
        self.ny = int(parts[3])
        self.nz = int(parts[4])
        self.x  = float(parts[5])
        self.y  = float(parts[6])
        self.z  = float(parts[7])
        self.dx = float(parts[8])
        self.dy = float(parts[9])
        self.dz = float(parts[10])
        return

    def dump(self):
        if self.fieldType == 'H':
            print('NH Near Magnetic Field:')
        else:
            print('NE Near Electric Field:')
        if self.coordType == 0:
            print('    Using Rectangular Coordinates')
            print('    Points desired in X:     %d' % self.nx)
            print('    Points desired in Y:     %d' % self.ny)
            print('    Points desired in Z:     %d' % self.nz)
            print('    First field point X:     %f m (%.2f in)' % (self.x, inches(self.x)))
            print('    First field point Y:     %f m (%.2f in)' % (self.y, inches(self.y)))
            print('    First field point Z:     %f m (%.2f in)' % (self.z, inches(self.z)))
            print('    Stepping increment in X  %f m (%.2f in)' % (self.dx, inches(self.dx)))
            print('    Stepping increment in Y  %f m (%.2f in)' % (self.dy, inches(self.dy)))
            print('    Stepping increment in Z  %f m (%.2f in)' % (self.dz, inches(self.dz)))
        elif self.coordType == 1:
            print('    Using Spherical Coordinates')
            print('    Points desired in r:         %d' % self.nx)
            print('    Points desired in phi:       %d' % self.ny)
            print('    Points desired in theta:     %d' % self.nz)
            print('    First field point r:         %f m (%.2f in)' % (self.x, inches(self.x)))
            print('    First field point phi:       %f degrees' % self.y)
            print('    First field point theta:     %f degrees' % self.z)
            print('    Stepping increment in r:     %f m (%.2f in)' % (self.dx, inches(self.dx)))
            print('    Stepping increment in phi:   %f degrees' % self.dy)
            print('    Stepping increment in theta: %f degrees' % self.dz)
        else:
            print('    Unexpected coordinate type argument %d' % self.coordType)
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
        self.opts = {}
        if verbose:
            self.opts['verbose'] = True

        self.comments = []
        self.geometry = Geometry()

        # Program control storage
        self.excitations      = []
        self.frequencies      = []
        self.ground           = None
        self.additionalGround = None
        self.loadings         = []
        self.networks         = []
        self.transmissionLines = []
        self.radiationPatterns = []
        self.nearFields       = []
        self.extendedKernel   = False
        self.interactionRange = None

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
            print('Comment: %s' % comm)
        return

    def CE(self, line):
        # End Comment
        sline = line.strip()
        comm = sline[3:]
        if comm != '':
            self.comments.append(comm)
        if 'verbose' in self.opts:
            print('Com End: %s' % comm)
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
        self.geometry.addWire(nw)
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
        self.geometry.end = ge
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
            print('GX Reflection in coordinate plane:')
            print('    Tag Increment:      %d' % TNI)
            if dstr[0] == '1':
                print('    Reflected along X axis')
            if dstr[1] == '1':
                print('    Reflected along Y axis')
            if dstr[2] == '1':
                print('    Reflected along Z axis')
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
            print('SP Surface Patch:')
            if NS == 0:
                # Arbitrary shape
                print('    Shape: Arbitrary')
            elif NS == 1:
                print('    Shape: Rectangular')
            elif NS == 2:
                print('    Shape: Triangular')
            elif NS == 3:
                print('    Shape: Quadrilateral')
            print('    X1     %f m (%.2f in)' % (X1, inches(X1)))
            print('    Y1     %f m (%.2f in)' % (Y1, inches(Y1)))
            print('    Z1     %f m (%.2f in)' % (Z1, inches(Z1)))
            print('    X2     %f m (%.2f in)' % (X2, inches(X2)))
            print('    Y2     %f m (%.2f in)' % (Y2, inches(Y2)))
            print('    Z2     %f m (%.2f in)' % (Z2, inches(Z2)))
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
            print('SM Multiple Patch Surface:')
            print('    Patches across X:     %d' % NX)
            print('    Patches across Y:     %d' % NY)
            print('    Corner 1 x:           %f m (%.2f in)' % (X1, inches(X1)))
            print('    Corner 1 y:           %f m (%.2f in)' % (Y1, inches(Y1)))
            print('    Corner 1 z:           %f m (%.2f in)' % (Z1, inches(Z1)))
            print('    Corner 2 x:           %f m (%.2f in)' % (X2, inches(X2)))
            print('    Corner 2 y:           %f m (%.2f in)' % (Y2, inches(Y2)))
            print('    Corner 2 z:           %f m (%.2f in)' % (Z2, inches(Z2)))
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
            print('SC Surface Multiple Patch Continuation:')
            print('    Corner 3 x:           %f m (%.2f in)' % (X3, inches(X3)))
            print('    Corner 3 y:           %f m (%.2f in)' % (Y3, inches(Y3)))
            print('    Corner 3 z:           %f m (%.2f in)' % (Z3, inches(Z3)))
        return

    #
    # Program Control
    #

    def CP(self, line):
        # Maximum Coupling Calculation
        if 'verbose' in self.opts:
            print('CP Maximum Coupling Calculation:')
            print('    Not supported in nec2c')
        return

    def EK(self, line):
        # Extended thin-wire kernel
        parts = line.split()
        I1 = int(parts[1])
        self.extendedKernel = (I1 == 0)
        if 'verbose' in self.opts:
            print('EK Extended thin-wire kernel:')
            if I1 == 0:
                print('    Initiate extended thin-wire kernel')
            elif I1 == -1:
                print('    Return to standard thin-wire kernel')
            else:
                print('    Unexpected value %d' % I1)
        return

    def EX(self, line):
        # Excitation
        ex = Excitation(line)
        self.excitations.append(ex)
        if 'verbose' in self.opts:
            ex.dump()
        return

    def FR(self, line):
        # Frequency
        self.frequencies.append(Frequency(line))
        if 'verbose' in self.opts:
            self.frequencies[-1].dump()
        return

    def GN(self, line):
        # Ground parameters
        self.ground = Ground(line)
        if 'verbose' in self.opts:
            self.ground.dump()
        return

    def KH(self, line):
        # Interaction Approximation Range
        parts  = line.split()
        self.interactionRange = int(parts[1])
        if 'verbose' in self.opts:
            print('KH Interaction Approximation Range:')
            print('    Approximation used for interactions over %d wavelengths' % self.interactionRange)
        return

    def LD(self, line):
        # Loading
        ld = Loading(line)
        if ld.type == -1:
            self.loadings = []
        else:
            self.loadings.append(ld)
        if 'verbose' in self.opts:
            ld.dump()
        return

    def GD(self, line):
        # Additional ground parameters
        self.additionalGround = AdditionalGround(line)
        if 'verbose' in self.opts:
            self.additionalGround.dump()
        return

    def NH(self, line):
        # Near Field NH
        self.nearFields.append(NearField(line, 'H'))
        if 'verbose' in self.opts:
            self.nearFields[-1].dump()
        return

    def NE(self, line):
        # Near Field NE
        self.nearFields.append(NearField(line, 'E'))
        if 'verbose' in self.opts:
            self.nearFields[-1].dump()
        return

    def NT(self, line):
        # Networks
        # Untested
        nt = Network(line)
        self.networks.append(nt)
        if 'verbose' in self.opts:
            nt.dump()
        return

    def NX(self, line):
        # Next Structure
        print('NX Next Structure:')
        #print('    %s' % line, end='')
        print('    Not supported in nec2c')
        return

    def RP(self, line):
        # Radiation Pattern
        self.radiationPatterns.append(RadiationPattern(line))
        if 'verbose' in self.opts:
            self.radiationPatterns[-1].dump()
        return

    def PQ(self, line):
        # Print control for charge on wires
        print('PQ Print control for charge on wires:')
        #print('    %s' % line, end='')
        print('    Not supported in nec2c')
        return

    def PT(self, line):
        # Page Title / Print Control for Current on Wires
        print('PT Page Title / Print Control for Current on Wires :')
        #print('    %s' % line, end='')
        print('    Not supported in nec2c')
        return

    def TL(self, line):
        # Transmission Line
        tl = TransmissionLine(line)
        self.transmissionLines.append(tl)
        if 'verbose' in self.opts:
            tl.dump()
        return

    def WG(self, line):
        # Write NGF file
        print('WG Write NGF file:')
        #print('    %s' % line, end='')
        print('    Not supported in nec2c')
        return

    def XQ(self, line):
        # Execute
        print('XQ Execute:')
        #print('    %s' % line, end='')
        print('    Ignored in nec2c')
        return

    def EN(self, line):
        # End of run
        if 'verbose' in self.opts:
            print('EN End of run')
        return

    def dump(self):
        for c in self.comments:
            print('Comment: %s' % c)
        self.geometry.dump()
        for fr in self.frequencies:
            fr.dump()
        for ex in self.excitations:
            ex.dump()
        for ld in self.loadings:
            ld.dump()
        if self.ground is not None:
            self.ground.dump()
        if self.additionalGround is not None:
            self.additionalGround.dump()
        if self.extendedKernel:
            print('EK Extended thin-wire kernel: Initiated')
        if self.interactionRange is not None:
            print('KH Interaction Approximation Range: %d wavelengths' % self.interactionRange)
        for nt in self.networks:
            nt.dump()
        for tl in self.transmissionLines:
            tl.dump()
        for nf in self.nearFields:
            nf.dump()
        for rp in self.radiationPatterns:
            rp.dump()
        print('EN End of run')
        return

def main():
    if len(sys.argv) < 2:
        print('Usage: %s filename [--verbose] [--dump]' % sys.argv[0])
        return

    verbose = '--verbose' in sys.argv or '-v' in sys.argv
    model = NecModel(verbose=verbose)

    inf = open(sys.argv[1])
    for line in inf.readlines():
        lparts = line.split()
        if not lparts:
            continue
        op = lparts[0].upper()
        if op in model.ops:
            model.ops[op](line)
        else:
            print('Unhandled -> ', line, end='')

    if '--dump' in sys.argv or '-d' in sys.argv:
        model.dump()

    return model

main()
