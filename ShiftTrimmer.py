import numpy as np

from pyraf import iraf,pyrafglobals
iraf.images()

### File constants
IMDIR = ''
IMROOT = 'HatP12b.Light.60S1X2.V.'
IMEXT = '' 

### CCD size
CCDx = 1092
CCDy = 736

def shiftTrim(files, xi, xf, yi, yf, dx=None, dy=None, test=False):
    """
    Given initial and final x and y values of shifted stars, will compute shift in x and y and 
    trim files to compensate for shifting image.

    Warning:
    - Assumes that movement across CCD is uniform and predictable in x and y (can be zero)
    - Discard bad images after trim is complete (requires complete series to trim accurately)
    - Only use reduced science images 
    - Pixels are discretely counted

    Parameters:
    ----------------------
    parameter: (dtype) [default (if optional)], information

    xi: (int), initial x value(s) for star(s)
    xf: (int), final x values(s) for star(s)
    yi: (int), initial y value(s) for star(s)
    yx: (int), final y value(s) for star(s)
    dx: (int) [None], overwrite x shift value (will ignore xi,xf)
    dy: (int) [None], overwrite y shift value (will ignore yi,yf)
    test: (boolean) [False], if True will only print pixel output and not trim files
    ----------------------
    """
    num = len(files) - 1
    
    if dx == None:
        dx = _dCalc(xi,xf,num,'x')
        
    if dy == None:
        dy = _dCalc(yi,yf,num,'y')
    
    print 'dx: %s, dy: %s' % (dx,dy)
    for i,f in enumerate(files):
        x = _pixelFinder(dx, i, abs(dx*num), 1, CCDx)
        y = _pixelFinder(dy, i, abs(dy*num), 1, CCDy)
        
        print "x: %s, y: %s, x*y: %s" % (x, y, ((x[1]-x[0])*(y[1]-y[0])))
        
        f_trim = f + '[%s:%s,%s:%s]' % (x[0],x[1],y[0],y[1])
        if test == False:
            iraf.imcopy(f_trim, f)
        

def fileFinder(minNum, maxNum, imdir=IMDIR, imroot=IMROOT, imext=IMEXT):
    """
    Finds files in range defined by minNum and maxNum using default class file roots and extension.

    Parameters:
    ----------------------
    parameter: (dtype) [default (if optional)], information

    minNum: (int), minimum file number
    maxNum: (int), maximum file number
    imdir: (string) [IMDIR], image directory
    imroot: (string) [IMROOT], image root (non-variable file name)
    imext: (string) [IMEXT], file extension
    ----------------------
    """
    batch = np.arange(minNum,maxNum+1,1)
    files = []

    for f in batch:
        files.append(imdir + imroot + str(f) + imext)

    return files

def _dCalc(i,f,num,axis):
    d = int((np.average(f - i) / num))
    
    if d == 0:
        if (f - i) > 0:
            d = 1
        elif (f - i) < 0:
            d = -1
        else:
            print 'No pixel shift detected in %s' % axis
            d = 0

    return d
    
def _pixelFinder(dp, i, dpmax, pmin, pmax):
    if dp > 0:
        return [pmin + dp*i, (pmax - dpmax) + dp*i]
    else: 
        return [dpmax + 1 + dp*i, pmax + dp*i]