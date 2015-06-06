import numpy as np

from pyraf import iraf,pyrafglobals
iraf.images()

### File constants
IMDIR = ''
IMROOT = 'HatP12b.Light.60S1X2.V.'
IMEXT = '' 

TRIMDIR = 'trimmedimages/'

### CCD size
CCDx = 1092
CCDy = 736

def shiftTrim(files, xi, xf, yi, yf, dx=None, dy=None, offset_x=0, offset_y=0, 
    fnew='', newdir=TRIMDIR, test=False):
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
    fnew: (string) [None], add string to new file name
    newdir: (string) [TRIMDIR], trimmed image directory
    test: (boolean) [False], if True will only print pixel output and not trim files
    ----------------------
    """
    xi = np.array(xi)
    xf = np.array(xf)
    yi = np.array(yi)
    yf = np.array(yf)
    
    num = len(files) - 1
    
    if dx == None:
        dx = _dCalc(xi,xf,num,'x')
        
    if dy == None:
        dy = _dCalc(yi,yf,num,'y')
    
    print 'dx: %s, dy: %s' % (dx,dy)
    for i,f in enumerate(files):

        x = _pixelFinder(dx, i, abs(dx*num), 1, CCDx, offset_x)
        y = _pixelFinder(dy, i, abs(dy*num), 1, CCDy, offset_y)
        
        print "x: %s, y: %s, x*y: %s" % (x, y, ((x[1]-x[0])*(y[1]-y[0])))
        
        f_trim = f + '[%s:%s,%s:%s]' % (x[0],x[1],y[0],y[1])

        if fnew != None:
            f = f + '.' + fnew

        if test == False:
            iraf.imcopy(f_trim, newdir + f)
        

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

def batchOperation(minNum, maxNum, batchsize, dx, dy, fnew='', offset_batchnum=0, newdir=TRIMDIR, imdir=IMDIR, 
    imroot=IMROOT, imext=IMEXT, test=False):

    mins, maxs, batchnum, remainder, batchnames = batchFinder(minNum, maxNum, batchsize, fnew, offset_batchnum)

    for i in range(batchnum):
        files = fileFinder(mins[i],maxs[i])

        if remainder != 0 and i == max(range(batchnum)):
            offset_x = _offsetCalc(dx, batchsize, remainder)
            offset_y = _offsetCalc(dy, batchsize, remainder)
            shiftTrim(files,0,0,0,0, dx=dx, dy=dy, offset_x=offset_x, offset_y=offset_y, fnew=batchnames[i], 
                newdir=newdir, test=test)
        else:
            shiftTrim(files,0,0,0,0, dx=dx, dy=dy, offset_x=0, offset_y=0, fnew=batchnames[i], newdir=newdir, test=test)

    return

def batchFinder(minNum, maxNum, batchsize, fnew, offset_batchnum=0):

    total = maxNum - minNum
    batchnum = (total + 1) / batchsize
    remainder = (total + 1) % batchsize

    mins = []
    maxs = []
    batchnames = []
    
    for i in range(batchnum):
        mins.append(minNum + batchsize*i)
        maxs.append(minNum + batchsize + batchsize*i - 1)
        batchnames.append(fnew + '.b' + str(offset_batchnum + i + 1))

    if remainder != 0:
        mins.append(max(maxs)+1)
        maxs.append(max(maxs)+ remainder)
        batchnum += 1
        batchnames.append(fnew + '.b' + str(offset_batchnum + batchnum))

    return mins, maxs, batchnum, remainder, batchnames

def _offsetCalc(d, batchsize, remainder):
    offset = batchsize - remainder
    return (offset*abs(d))

def _dCalc(i,f,num,axis):
    d = int((np.average(f - i) / num))
    
    if d == 0:
        if np.average(f - i) > 0:
            d = 1
        elif np.average(f - i) < 0:
            d = -1
        else:
            print 'No pixel shift detected in %s' % axis
            d = 0

    return d
    
def _pixelFinder(dp, i, dpmax, pmin, pmax, offset):
    if dp > 0:
        return [offset + pmin + dp*i, (pmax - dpmax) + dp*i]
    else: 
        return [dpmax + 1 + dp*i, pmax - offset + dp*i]