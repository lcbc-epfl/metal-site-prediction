import time
import multiprocessing
from multiprocessing import Pool

import torch

from moleculekit.molecule import Molecule
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors
from moleculekit.tools.atomtyper import prepareProteinForAtomtyping
from moleculekit.tools.preparation import systemPrepare


class AtomtypingError(Exception):
    pass

class StructureCleaningError(Exception):
    pass

class ProteinPrepareError(Exception):
    pass

class VoxelizationError(Exception):
    pass

metal_atypes = ('MG', 'ZN', 'MN', 'CA', 'FE', 'HG', 'CD', 'NI', 'CO', 'CU', 'K', 'LI', 'Mg', 'Zn', 'Mn', 'Ca', 'Fe', 'Hg', 'Cd', 'Ni', 'Co', 'Cu', 'Li')

def voxelize_single_notcentered(env):
    """voxelize 1 structure, executed on a single CPU 
    Using 7 of the 8 channels supplied by moleculekit(excluding metals)
    Additionally it uses all the metalbinding residues as channel

    Parameters
    ----------
    env : tuple
        Tuple of the form (prot, idx)

    Returns
    -------
    voxels : torch.tensor
        Voxelized structure with 8 channels (8,20,20,20)
    prot_centers : list
        List of the centers of the voxels (20x20x20,3)
    prot_n : list
        List of the number of voxels in each voxel (20x20x20)
    prot : moleculekit.Molecule
        Moleculekit molecule
    """
    prot, id = env

    c = prot.get('coords',sel=f"resid {id} and name CA")

    size=[20,20,20] # size of box
    voxels = torch.zeros(8,20,20,20)

    try:
        prot_vox, prot_centers, prot_N = getVoxelDescriptors(prot, buffer=1, center=c, boxsize=size,method="C", validitychecks=False)
    except:
        return None
    nchannels = prot_vox.shape[1]
    prot_vox_t = prot_vox.transpose().reshape([1, nchannels, prot_N[0], prot_N[1], prot_N[2]]).copy()
    metalcoordination = prot.atomselect('(resname HID HIS HIE HIP ASP ASH GLU GLH CYS and sidechain) or (protein and name O N)')
    metalcoordination = metalcoordination.reshape(metalcoordination.shape[0],1)
    coord_vox, _, _ = getVoxelDescriptors(prot, center=c, voxelsize=1, userchannels=metalcoordination, boxsize=[20,20,20], method="C", validitychecks=False)
    coord_vox_t = coord_vox.transpose().reshape([1, prot_N[0], prot_N[1], prot_N[2]])
    
    prot_vox_t = torch.from_numpy(prot_vox_t)
    indices = torch.LongTensor([0,1,2,3,4,5,7])
    voxels[0:7]=torch.index_select(prot_vox_t,1, indices)

    voxels[7]=torch.from_numpy(coord_vox_t)
    
    return (voxels, prot_centers, prot_N, prot.copy())


def processStructures(pdb_file, resids,clean=True):
    """Process a pdb file and return a list of voxelized boxes centered on the residues
    
    Parameters
    ----------
    pdb_file : str
        Path to pdb file
    resids : list
        List of resids to center the voxels on
    clean : bool
        If True, remove all non-protein residues from the pdb file
    
    Returns
    -------
    voxels : torch.Tensor
        Voxelized boxes with 8 channels (N, 8,20,20,20)
    prot_centers_list : list
        List of the centers of the voxels (N*20**20*20,3)
    prot_n_list : list
        List of the number of voxels in each box (N,3)
    envs: list
        List of tuples (prot, idx) (N)
    """
    
    start_time_processing = time.time()

    # load molecule using MoleculeKit
    try:
        prot = Molecule(pdb_file)
    except:
        raise IOError('could not read pdbfile')

    if clean:
        prot.filter('protein')
        prot.renumberResidues()
        
        # prepare for atomtyping
        try:
            # prot= proteinPrepare(prot) switch to new systemprepare
            prot= systemPrepare(prot, titration=False, verbose=False)
            prot.filter("protein or water or element {}".format(" ".join(metal_atypes))) #prot.remove("name OXT H2 H3") # this is a necessary hack because otherwise atomtyping might fail
        except: 
            raise StructureCleaningError('structure could not be cleaned')

    try:
        # prot = prepareProteinForAtomtyping(prot)
        prot = prepareProteinForAtomtyping(prot, protonate=False, verbose=False)
        prot.filter("protein or water or element {}".format(" ".join(metal_atypes)))
    except: 
        raise ProteinPrepareError('failed to prepare protein for atomtyping')

   
    environments = []
    for idx in resids:
        try:
            environments.append((prot.copy(), idx))
        except:
            print('ignoring '+idx)

    

    prot_centers_list = []
    prot_n_list = []
    envs = []

    cpuCount=multiprocessing.cpu_count()
    p = Pool(cpuCount)

    results = p.map(voxelize_single_notcentered, environments)  # prot_centers, prot_N

    # remove None results
    results = [x for x in results if x is not None]

    voxels=torch.empty(len(results),8,20,20,20)

    if len(results)==0:
        raise VoxelizationError('something went wrong with the voxelization, check that there are no discontinuities in the protein')

    vox_env, prot_centers_list, prot_n_list, envs = zip(*results) 

    for i,vox_env in enumerate(vox_env):
        voxels[i]=vox_env

    print(f"Voxelization took  {time.time() - start_time_processing:.3f} seconds ")

    return voxels, prot_centers_list, prot_n_list, envs