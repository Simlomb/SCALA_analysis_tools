#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
## Filename:          plan_scala.py
## Version:           $Revision: 1.0 $
## Description:
## Author:            Nicolas Chotard <nchotard@ipnl.in2p3.fr>
## Author:            $Author: nchotard $
## Created at:        $Date: 03-11-0015 14:59:57 $
## Modified at:       16-11-2015 10:59:22
## $Id: plan_scala.py, v 1.0, 03-11-0015 14:59:57 nchotard Exp $
################################################################################

"""
"""

__author__ = "Nicolas Chotard <nchotard@ipnl.in2p3.fr>"
__version__ = '$Id: plan_scala.py, v 1.0, 03-11-0015 14:59:57 chotard Exp $'

import os
import glob
import optparse

import numpy as N
from processing.process.models import Process

CODE_NAME = os.path.basename(__file__)

def buildName(name, pre="", ext=""):
    """Build filename from name: remove path and extension, prepend prefix and
    append extension)."""

    p, f = os.path.split(name)          # Remove path
    f = os.path.splitext(f)[0]          # Remove extension
    if ext:                             # Add an extension
        f = os.path.extsep.join((f, ext))
    if pre:
        f = pre + f                     # Prepend prefix
    return f
 
class ScalaCalib:

    def __init__(self, wcube, ffcube):
        self.channel = 'B' if '_B' in wcube else 'R'
        if self.channel == 'B':
            self. wrange = "3300,5150" # Don't keep crap at ends!
        else:
            self.wrange = "5100,10150" # Extended range for updated snifsmasks
        self.wcube = wcube # input cube
        self.ffcube = ffcube # flat
        self.ccp_cube = 'C'+self.wcube[1:] # cosmic cleaned (and flatfielded)
        self.uf_cube = 'UF'+self.wcube[1:] # cosmic cleaned (un-flatfielded)
        self.tuf_cube = 'TUF'+self.wcube[1:] # truncated version
        self.e3d_cube = "e3d_" + buildName(self.uf_cube) # Euro 3d cube
        self.threed_cube = "3d_" + buildName(self.uf_cube) # 3d cube
        self.e3d_tcube = "e3d_" + buildName(self.tuf_cube) # Euro 3d cube
        self.threed_tcube = "3d_" + buildName(self.tuf_cube) # 3d cube

    def quickCalib(self, verbose=True):
        """Build the `quick_calib` command, with appropriate flat, flux and
        domain options."""
        cmd = "quick_calib -i %s" % (buildName(self.wcube)) #incube
        cmd += ' -l ' + self.ffcube # flat
        cmd += ' -x None' # no flux calibration
        if verbose:
            print "# Run quick_calib to get the flatfieled/cosmic-cleaned cube"
            print cmd
        return cmd

    def unflatfield(self, verbose=True):
        cmd1 = "wr_desc -file %s -desc SFCLASS -val 1 -type int -quiet" % \
              self.ccp_cube
        if verbose:
            print "# Set SFCLASS to disable Fclass checks"
            print cmd1
        cmd2 = "apply_lfff -in %s -out %s -lfff %s -noask -reverse" % \
              (self.ccp_cube, self.uf_cube, self.ffcube)
        if verbose:
            print "# Un-flatfield the flatfieled/cosmic-cleaned cube"
            print cmd2
        return cmd1, cmd2

    def truncate(self, verbose=True):
        cmd = "truncate_cube -in %s -out %s -wave %s" % \
              (self.uf_cube, self.tuf_cube, self.wrange)
        if verbose:
            print "# Truncate the cosmic-cleaned and un-flatfieled cube"
            print cmd
        return cmd

    def convertfiles(self, verbose=True):
        """Build the `convert_file` command to produce Euro3D and 3D cubes."""

        cmd1 = 'convert_file -in %s -out %s ' \
              '-inputformat "tiger+fits" -outputformat "euro3d" -quiet' % \
              (self.uf_cube, self.e3d_cube)
        if verbose:
            print "# Convert the untruncated cube from tig to e3d"
            print cmd1
        cmd2 = 'convert_file -in %s -out %s ' \
              '-inputformat "tiger+fits" -outputformat "euro3d" -quiet' % \
              (self.tuf_cube, self.e3d_tcube)
        if verbose:
            print "# Convert the truncated cube from tig to e3d"
            print cmd2

        cmd3 = "e3dto3d.py -o %s.fits %s.fits" % (self.threed_cube, \
                                                  self.e3d_cube)
        if verbose:
            print "# Convert the untruncated cube from e3d to 3d"
            print cmd3
        cmd4 = "e3dto3d.py -o %s.fits %s.fits" % (self.threed_tcube, \
                                                  self.e3d_tcube)
        if verbose:
            print "# Convert the truncated cube from e3d to 3d"
            print cmd4
        return cmd1, cmd2, cmd3, cmd4

    def write_commands(self):
        self.quickCalib()
        self.unflatfield()
        self.truncate()
        self.convertfiles()

    def run_commands(self):
        cmds = []
        cmds.append([self.quickCalib(verbose=False)])
        cmds.append(self.unflatfield(verbose=False))
        cmds.append([self.truncate(verbose=False)])
        cmds.append(self.convertfiles(verbose=False))
        for cmd in N.concatenate(cmds):
            print "\n" + cmd
            print os.system(cmd)

def get_scala_cubes(night, version=202, status=2, fclass=65, run=None):
    """
    Get all the SCALA data from a given night
    """
    run = '' if run is None else str(run).zfill(3)
    scala_cube = Process.objects.filter(Version=version, Fclass=fclass,
                                        IdProcess__startswith=night+run)
    return scala_cube

def copy_scala_cubes(cubes):
    """
    Copy all the scala cubes from the DB to the local directory
    """
    frames = []
    for cube in cubes:
        fcube = cube.Parent_FK.get(Fclass=300).AuxFiles_FK.get(XFclass=3)
        wcube = cube.AuxFiles_FK.get(XFclass=2)
        wtable = cube.AuxFiles_FK.get(XFclass=1)
        os.system("cp -v %s ./%s" % (wcube.FullPath, wcube.OName.lstrip('D')))
        os.system("cp -v %s ./%s" % (wtable.FullPath, wtable.OName.lstrip('D')))
        os.system("cp -v %s ./%s" % (fcube.FullPath, fcube.OName.lstrip('D')))
        frames.append(ScalaCalib(wcube.OName.lstrip('D'),
                                 fcube.OName.lstrip('D')))
    return frames

def get_clap_data(night, run=None):
    """
    Has to be changed when the CLAP data will be stored into the DB. They are
    for now in /sps/snovae/SRBregister/SCALA/SCALA_Commissioning/CLAP_data
    """
    clapdir = '/sps/snovae/SRBregister/SCALA/SCALA_Commissioning/CLAP_data'
    if run is None:
        nightdir = clapdir+'/%2d/%3d/*.fits'%(int(night[:2]), int(night[2:]))
    else:
        nightdir = clapdir+'/%2d/%3d/*_%s_*.fits'%(int(night[:2]),
                                                   int(night[2:]),
                                                   str(run).zfill(3))
    clapdata = glob.glob(nightdir)
    return clapdata

def copy_clap_data(claps):
    """
    Copy a list of claps data to the local directory
    """
    for clap in claps:
        os.system("cp -v %s ." % clap)

def run_preprocessing(frames):
    """Run all the preprocessing steps on the w-caloibrated cubes"""
    for frame in frames:
        frame.run_commands()

## TESTS ##
def get_test_data(idprocess='151511170082'):
    """
    Get the SCALA test data from the DB
    """
    scala_cube = Process.objects.get(Version=202, Fclass=65,
                                     IdProcess__startswith=idprocess)
    flatfield_cube = scala_cube.Parent_FK.get(Fclass=300).AuxFiles_FK.get(XFclass=3)
    wcalib_cube = scala_cube.AuxFiles_FK.get(XFclass=2)
    wcalib_table = scala_cube.AuxFiles_FK.get(XFclass=1)
    os.system("cp -v %s ./%s" % (wcalib_cube.FullPath,
                                 wcalib_cube.OName.lstrip('D')))
    os.system("cp -v %s ./%s" % (wcalib_table.FullPath,
                                 wcalib_table.OName.lstrip('D')))
    os.system("cp -v %s ./%s" % (flatfield_cube.FullPath,
                                 flatfield_cube.OName.lstrip('D')))
    return wcalib_cube.OName.lstrip('D'), flatfield_cube.OName.lstrip('D')

def test():
    """
    Test the code
    """
    wcube, ffcube = get_test_data()
    cube = ScalaCalib(wcube, ffcube)
    return cube

def get_calib_files(year=15, day=159, channel=2):
    """
    Return telluric correction, extinction and flux solution processes for
    a given night (year, day) and channel
    """
    FclassTell = 620  # Telluric correction (w/ XFclassP)
    FclassExt = 625   # Extinction fclass (w/ XFclassS|XFclassM)
    FclassFxSol = 630 # Flux solution fclass (w/ XFclassS|XFclassM)
    XFclassP = 600    # XFclass multi-std photometric case
    # Nightly multi-std telluric correction from PMS
    tc = Process.objects.filter(Pose_FK__Exp_FK__Run_FK__Year=year,
                                Pose_FK__Exp_FK__Run_FK__Day=day,
                                Fclass=FclassTell, XFclass=XFclassP,
                                Version=202)
    # Photometric nightly multi-std extinction, from PMS
    ex = Process.objects.filter(Pose_FK__Exp_FK__Run_FK__Year=year,
                                Pose_FK__Exp_FK__Run_FK__Day=day,
                                Fclass=FclassExt, XFclass=XFclassP,
                                Version=202)
    # Look for flux solution for current night
    fs = Process.objects.filter(Pose_FK__Exp_FK__Run_FK__Year=year,
                                Pose_FK__Exp_FK__Run_FK__Day=day,
                                Fclass=FclassFxSol, XFclass=XFclassP,
                                Channel=channel, Version=202)
    return tc, ex, fs

# MAIN #########################################################################

if __name__ == "__main__":

    # Options ==================================================================
    usage = " usage: [%prog] [options] -n night"

    parser = optparse.OptionParser(usage, version=__version__)
    parser.add_option('-n', '--night', type='string', help="Night to work on")
    parser.add_option('-r', '--run', type='string', help="Run to work on")
    opts, args = parser.parse_args()

    # make sure the night has the right format
    night = str(opts.night).replace('_', '')

    # get the data
    scala_cubes = get_scala_cubes(opts.night, run=opts.run)
    clap_data = get_clap_data(opts.night, run=opts.run)

    # make local copy of the data
    frames = copy_scala_cubes(scala_cubes)
    copy_clap_data(clap_data)

    # run the preprocessing
    run_preprocessing(frames)

    # Hard coded path: BAD, but working for now
    path = '/afs/in2p3.fr/group/snovae/snprod/SnfProd/nchotard/plan_scala'

    # throuput estimate
    cmd = "python %s/SCALA_analysis_tools/SCALA_main.py " % path
    cmd += "-s %s -c %s" % (",".join(glob.glob('3d_U*.fits')),
                            ",".join(glob.glob('SC*.fits')))
    print cmd
    os.system("cp %s/SCALA_analysis_tools/*.txt ." % path)
    os.system(cmd)

    # flux calibration of scala cubes
    cmd = "python %s/SCALA_analysis_tools/SCALA_main.py " % path
    cmd += "-s %s -c %s -t False" % (",".join(glob.glob('3d_U*.fits')),
                                     ",".join(glob.glob('SC*.fits')))
    print cmd
    os.system("cp %s/SCALA_analysis_tools/*.txt .")
    os.system(cmd)

# End of plan_scala.py
