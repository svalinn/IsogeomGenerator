from pyne import variancereduction as vr
from pyne.mcnp import Meshtal, Wwinp

if __name__ == "__main__":
    
    f = "../run/meshtal"
    M = Meshtal(f)
    tally_number = M.tally.keys()[0]
    meshtally = M.tally[tally_number]
    tags = meshtally.tag_names
    
    res_tag = tags[0]
    err_tag = tags[1]
    
    vr.magic(meshtally, res_tag, err_tag)

    meshtally.write_hdf5("mesh-ww-tags.h5m")
    
    wwinp = Wwinp()
    wwinp.read_mesh(meshtally.mesh)
    
    # need this if you want to run MCNP again
    wwinp.write_wwinp("wwinp")
