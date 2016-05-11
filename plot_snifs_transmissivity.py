import numpy                         as N
import matplotlib.pyplot             as P
import pickle                        as pk
import optparse

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file",
                      help="txt file to open and plot",
                      )
    opts,args = parser.parse_args()
    cal= N.loadtxt(str(opts.file))
    
    P.plot(cal[:,0], cal[:,1], '-ob', label='transmissivity_%s' %str(opts.file)[-8])
    P.legend(loc='best', fancybox=True)
    P.xlabel("Wavelength ($\AA$)")
    P.ylabel("Transmissivity")
    P.savefig('transmissivity_%s.png' %str(opts.file)[:-7], dpi=300)
    
