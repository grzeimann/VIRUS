
import cosmics
import numpy as np
import os.path

__author__ = 'Dustin Davis'

def remove_cosmics(sci_file, prefix='c'):
    path = os.path.dirname(sci_file)
    if len(path) > 0:
        path += "/"
    base = os.path.basename(sci_file)

    new_sci_file = path + prefix + base
    new_err_file = path + 'e.' + prefix + base

    array, header = cosmics.fromfits(sci_file)
    e_array, e_header = cosmics.fromfits(path + 'e.' + base)

    # split in half so can use individual readnoise
    array_l, array_u = np.split(array, 2, 1)

    readnoise_u = 1.0
    readnoise_l = 1.0
    try: # average read noise
        readnoise_u = float(header['RDNOISEU'])
        readnoise_l = float(header['RDNOISEL'])
    except:
        pass

    # Build the object :
    # There are other options, check the manual...
    # sigclip --- want sharp edge; noise is around 10, sources pretty much never above 120-140, so 20x on noise is safe
    # sigfrac -- pixels neighboring cosmics ... because sigclip is high and cosmics can be hot (10K+) this needs to be small
    # objlim -- needs to be small, real sources very faint vs cosmics (this is a contrast between cosmics and source)
    cu = cosmics.cosmicsimage(array_u, gain=1.0, readnoise=readnoise_u, sigclip=25.0, sigfrac=0.001, objlim=0.001, satlevel=-1.0)
    cl = cosmics.cosmicsimage(array_l, gain=1.0, readnoise=readnoise_l, sigclip=25.0, sigfrac=0.001, objlim=0.001, satlevel=-1.0)

    # Run the removal
    # note: sometimes get warnings about invalid values encountered. This appears due to NaN in the fits error files (noise)
    cu.run(maxiter=4)
    cl.run(maxiter=4)

    # recombine arrays
    # use the mask to make our edits
    c = np.where(cl.mask == True)
    for x, y in zip(c[0], c[1]):
        array[x][y] = 0.0
        e_array[x][y] = -1.0

    x_off,y_off = array_u.shape
    c = np.where(cu.mask == True)
    for x, y in zip(c[0], c[1]):
        array[x][y+y_off] = 0.0
        e_array[x][y+y_off] = -1.0

    # Write the cleaned image into a new FITS file, conserving the original header :
    header['HISTORY']='removed cosmics (LACosmics)'
    e_header['HISTORY'] = 'removed cosmics (LACosmics)'

    cosmics.tofits(new_sci_file, array, header)
    cosmics.tofits(new_err_file, e_array, e_header)

    return
