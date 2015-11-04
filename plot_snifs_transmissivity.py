import numpy                         as N
import matplotlib.pyplot             as P
import pickle                        as pk
import optparse

if __name__ == '__main__':
    parser = optparse.OptionParser()
    parser.add_option("-f", "--file",
                      help="pickle file to open and plot",
                      )
    opts,args = parser.parse_args()
    with open(str(opts.file), 'rb') as transm_info:
        transmittivity = pk.load(transm_info)
    new_trans = N.zeros((len(transmittivity[0,:,0,0]),15,15))
    new_lambd = N.zeros((len(transmittivity[0,:,0,0]),15,15))
    for i in range(15):
        for j in range(15):
            new_trans[:,i,j] = transmittivity[1,:,i,j]
            new_lambd[:,i,j] = transmittivity[0,:,i,j]


    transm = N.mean(new_trans, axis=1)
    transm = N.mean(transm, axis=1)
    lbd = N.mean(new_lambd, axis=1)
    lbd = N.mean(lbd, axis=1)
    P.plot(lbd, transm, '-ob', label='transmissivity_%s' %str(opts.file)[-8])
    P.legend(loc='best', fancybox=True)
    P.xlabel("Wavelength ($\AA$)")
    P.ylabel("Transmissivity")
    P.savefig('transmissivity_%s.png' %str(opts.file)[:-7], dpi=300)
    
