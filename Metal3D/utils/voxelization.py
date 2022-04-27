import time
import multiprocessing
from multiprocessing import Pool

import torch
import numpy as np

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


metal_atypes = (
    "MG",
    "ZN",
    "MN",
    "CA",
    "FE",
    "HG",
    "CD",
    "NI",
    "CO",
    "CU",
    "K",
    "LI",
    "Mg",
    "Zn",
    "Mn",
    "Ca",
    "Fe",
    "Hg",
    "Cd",
    "Ni",
    "Co",
    "Cu",
    "Li",
)


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

    c = prot.get("coords", sel=f"index {id} and name CA")

    size = [16, 16, 16]  # size of box
    voxels = torch.zeros(8, 32, 32, 32)

    try:
        hydrophobic = prot.atomselect("element C")
        hydrophobic = hydrophobic.reshape(hydrophobic.shape[0], 1)

        aromatic = prot.atomselect(
            "resname HIS HIE HIP HID TRP TYR PHE and sidechain and not name CB and not hydrogen"
        )
        aromatic = aromatic.reshape(aromatic.shape[0], 1)

        metalcoordination = prot.atomselect(
            "(name ND1 NE2 SG OE1 OE2 OD2) or (protein and name O N)"
        )
        metalcoordination = metalcoordination.reshape(metalcoordination.shape[0], 1)

        hbondacceptor = prot.atomselect(
            "(resname ASP GLU HIS HIE HIP HID SER THR MSE CYS MET and name ND2 NE2 OE1 OE2 OD1 OD2 OG OG1 SE SG) or name O"
        )
        hbondacceptor = hbondacceptor.reshape(metalcoordination.shape[0], 1)

        hbonddonor = prot.atomselect(
            "(resname ASN GLN ASH GLH TRP MSE SER THR MET CYS and name ND2 NE2 NE1 SG SE OG OG1) or name N"
        )
        hbonddonor = hbonddonor.reshape(metalcoordination.shape[0], 1)

        positive = prot.atomselect(
            "resname LYS ARG HIS HIE HIP HID and name NZ NH1 NH2 ND1 NE2 NE"
        )
        positive = positive.reshape(positive.shape[0], 1)

        negative = prot.atomselect("(resname ASP GLU ASH GLH and name OD1 OD2 OE1 OE2)")
        negative = negative.reshape(negative.shape[0], 1)

        occupancy = prot.atomselect("protein and not hydrogen")
        occupancy = occupancy.reshape(occupancy.shape[0], 1)
        userchannels = np.hstack(
            [
                hydrophobic,
                aromatic,
                metalcoordination,
                hbondacceptor,
                hbonddonor,
                positive,
                negative,
                occupancy,
            ]
        )
        prot_vox, prot_centers, prot_N = getVoxelDescriptors(
            prot,
            center=c,
            userchannels=userchannels,
            boxsize=size,
            voxelsize=0.5,
            validitychecks=False,
        )
    except:
        raise VoxelizationError(f"voxelization of {id} failed")
    nchannels = prot_vox.shape[1]
    prot_vox_t = (
        prot_vox.transpose()
        .reshape([1, nchannels, prot_N[0], prot_N[1], prot_N[2]])
        .copy()
    )

    voxels = torch.from_numpy(prot_vox_t)
    return (voxels, prot_centers, prot_N, prot.copy())


def processStructures(pdb_file, resids, clean=True):
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
        Voxelized boxes with 8 channels (N, 8,32,32,32)
    prot_centers_list : list
        List of the centers of the voxels (N*32**32*32,3)
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
        raise IOError("could not read pdbfile")

    if clean:
        prot.filter("protein and not hydrogen")

    environments = []
    for idx in resids:
        try:
            environments.append((prot.copy(), idx))
        except:
            print("ignoring " + idx)

    prot_centers_list = []
    prot_n_list = []
    envs = []

    results = [voxelize_single_notcentered(x) for x in environments]

    voxels = torch.empty(len(results), 8, 32, 32, 32, device="cuda")

    vox_env, prot_centers_list, prot_n_list, envs = zip(*results)

    for i, vox_env in enumerate(vox_env):
        voxels[i] = vox_env

    print(f"Voxelization took  {time.time() - start_time_processing:.3f} seconds ")

    return voxels, prot_centers_list, prot_n_list, envs
