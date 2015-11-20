#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
## Filename:          plan_scala.py
## Version:           $Revision: 1.0 $
## Description:       Run the SCALA analysis on a given night
## Author:            Nicolas Chotard <nchotard@ipnl.in2p3.fr>
## Author:            $Author: nchotard $
## Created at:        $Date: 03-11-0015 14:59:57 $
## Modified at:       20-11-2015 11:45:57
## $Id: plan_scala.py, v 1.0, 03-11-0015 14:59:57 nchotard Exp $
################################################################################

"""
Run the SCALA analysis for a given night (run)
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
    # Remove path
    p, f = os.path.split(name)
    # Remove extension
    f = os.path.splitext(f)[0]
    if ext: # Add an extension          
        f = os.path.extsep.join((f, ext))
    if pre: # Prepend prefix
        f = pre + f                     
    return f

class ScalaCube:

    def __init__(self, scube):
        """
        Initialize a scala cube from a scala process (Fclass=65, XFclass=0)
        with parents and auxiliary files
        """
        self.scala_cube = scube
        self.channel = scube.Channel
        self.flat = scube.Parent_FK.get(Fclass=300).AuxFiles_FK.get(XFclass=3)
        self.wavelength_calibrated_table = scube.AuxFiles_FK.get(XFclass=1)
        self.wavelength_calbrated_cube = scube.AuxFiles_FK.get(XFclass=2)
        self.reduced_cube = scube.AuxFiles_FK.get(XFclass=3)
        
    def copy_files(self, path='./'):
        """
        Copu all the needed files to a given path
        """
        if not path.endswith('/'):
            path += '/'
        # the flat (tig)
        os.system("cp -v %s%s %s" % (path,self.flat.File_FK.FullPath,
                                     self.flat.File_FK.OName.lstrip('D')))
        # the wavelength calibrated table (tig)
        os.system("cp -v %s%s %s" % (path, self.flat.File_FK.FullPath,
                                     self.flat.File_FK.OName.lstrip('D')))
        # the wavelength calibrated cueb (tig)
        os.system("cp -v %s%s %s" % (path, self.flat.File_FK.FullPath,
                                   self.flat.File_FK.OName.lstrip('D')))
        # the reduced cube (tig)
        os.system("cp -v %s%s %s" % (path, self.flat.File_FK.FullPath,
                                     self.flat.File_FK.OName.lstrip('D')))
 
class ScalaCalib:

    def __init__(self, wcube, ffcube):
        self.channel = 'B' if '_B' in wcube else 'R'
        if self.channel == 'B':
            # Don't keep crap at ends!
            self. wrange = "3300,5150" 
        else:
            # Extended range for updated snifsmasks
            self.wrange = "5100,10150"
        # input cube
        self.wcube = wcube
        # flat
        self.ffcube = ffcube
        # cosmic cleaned (and flatfielded)
        self.ccp_cube = 'C'+self.wcube[1:]
        # cosmic cleaned (un-flatfielded)
        self.uf_cube = 'UF'+self.wcube[1:]
        # truncated version
        self.tuf_cube = 'TUF'+self.wcube[1:]
        # Euro 3d cube
        self.e3d_cube = "e3d_" + buildName(self.uf_cube)
        # 3d cube
        self.threed_cube = "3d_" + buildName(self.uf_cube)
        # Euro 3d cube
        self.e3d_tcube = "e3d_" + buildName(self.tuf_cube)
        # 3d cube
        self.threed_tcube = "3d_" + buildName(self.tuf_cube)

    def quickCalib(self, verbose=True):
        """
        Build the `quick_calib` command, with appropriate flat, flux and domain
        options.
        """
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

class ScalaFluxCalib:

    def __init__(self, scube):
        """
        Initialize the ScalaFluxCalib class with a SCALA cube (F/XFclass=65/0)
        """
        self.scube = scube

    def truncate_cube(self):
        print "# Truncate cube"
        cmd = "truncate_cube -in %s -out %s -wave 5100,9700 -quiet -noask" % ('TCP15_159_037_007_60_R.tig',
                                                                              'TTCP15_159_037_007_60_R.tig')

    def change_header(self):
        print "# Change airmass in the header"
        cmd1 = "wr_desc -file %s -desc AIRMASS -type double -val 0" % 'TTCP15_159_037_007_60_R.tig'
        cmd2 = "wr_desc -file %s -desc IAIRMASS -type double -val `rd_desc -file %s -desc AIRMASS -quiet`" % \
               ('TTCP15_159_037_007_60_R.tig', 'TCP15_159_037_007_60_R.tig')

    def flux_calibrate(self):
        print "# Apply flux calibration"
        cmd = "apply_flux -flux %s " % "fxSolP_15_159_R.fits"
        cmd += "-extinct %s,LAMBDA,EXT " % "extP_15_159.fits"
        cmd += "-in %s" % "TTCP15_159_037_007_60_R.tig "
        cmd += "-out %s" % "XTCP15_159_037_007_60_R.tig"

    def clean_table(self):
        print "Clean the table"
        cmd = "sel_table -in %s -sel all -quiet" % "CP15_159_037_007_60_R.fits"

    def convert_files(self):
        print "# Convert to e3d"
        cmd = "convert_file -inputformat 'tiger+fits' -outputformat 'euro3d' -quiet "
        cmd += "-in %s " % "XTCP15_159_037_007_60_R.tig"
        cmd += "-out %s" % "e3d_XTCP15_159_037_007_60_R.fits"
        
        print "# Convert to 3d"
        cmd = "e3dto3d.py -o %s %s" % ("3d_XTCP15_159_037_007_60_R.fits",
                                     "e3d_XTCP15_159_037_007_60_R.fits")

        
        #echo "# Copy files"
        #cp -v /sps/snovae/SRBregister/Prod/02-02/15/159/15_159_014_003_2_625_600_02-02_000.fits extP_15_159.fits
        #cp -v /sps/snovae/SRBregister/Prod/02-02/15/159/15_159_014_003_2_630_600_02-02_000.fits fxSolP_15_159_R.fits
        #cp -v /sps/snovae/SRBregister/Prod/02-02/15/159/15_159_037_007_2_065_003_02-02_000.tig TCP15_159_037_007_60_R.tig
        #cp -v /sps/snovae/SRBregister/Prod/02-02/15/159/15_159_037_007_2_065_001_02-02_000.fits CP15_159_037_007_60_R.fits

    def write_commands(self):
        self.truncate_cube()
        self.change_header()
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

def run_flux_calibration(frames, ex_B, fs_B, ex_R, fs_R):
    pass

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
    Return extinction and flux solution processes for
    a given night (year, day) and channel
    """
    FclassExt = 625   # Extinction fclass (w/ XFclassS|XFclassM)
    FclassFxSol = 630 # Flux solution fclass (w/ XFclassS|XFclassM)
    XFclassP = 600    # XFclass multi-std photometric case
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
    return ex, fs

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
    scala_cubes = get_scala_cubes(night, run=opts.run)
    clap_data = get_clap_data(night, run=opts.run)

    # make local copy of the data
    #calib_frames, flux_calib_frame = copy_scala_cubes(scala_cubes)
    calib_frames = copy_scala_cubes(scala_cubes)
    copy_clap_data(clap_data)

    # run the preprocessing
    #run_preprocessing(scala_cubes)
    run_preprocessing(calib_frames)

    """
    # get the needed files for flux calibration
    ex_B, fs_B = get_calib_files(year=15, day=159, channel=4)
    ex_R, fs_R = get_calib_files(year=15, day=159, channel=2)

    # run the flux calibration
    run_flux_calibration(frames)
    """

    # Hard coded path: BAD, but working for now
    path = '/afs/in2p3.fr/group/snovae/snprod/SnfProd/nchotard/plan_scala'

    # throuput estimate
    cmd = "python %s/SCALA_analysis_tools/SCALA_main.py " % path
    cmd += "-s %s -c %s" % (",".join(glob.glob('3d_U*.fits')),
                            ",".join(glob.glob('SC*.fits')))
    print cmd
    os.system("cp -f %s/SCALA_analysis_tools/*.txt ." % path)
    os.system(cmd)

    """
    # flux calibration of scala cubes
    cmd = "python %s/SCALA_analysis_tools/SCALA_main.py " % path
    cmd += "-s %s -c %s -t False" % (",".join(glob.glob('3d_U*.fits')),
                                     ",".join(glob.glob('SC*.fits')))
    print cmd
    os.system("cp %s/SCALA_analysis_tools/*.txt .")
    os.system(cmd)
    """

# End of plan_scala.py
