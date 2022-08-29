import os
import multiprocessing
from multiprocessing import Pool
from turtle import width

import py3Dmol
import numpy as np

from moleculekit.molecule import Molecule
from scipy.spatial import KDTree
from sklearn.cluster import AgglomerativeClustering
from tqdm.contrib.concurrent import process_map


def create_grid_fromBB(boundingBox, voxelSize=1):
    """Create a grid from a bounding box.

    Parameters
    ----------
    boundingBox : list
        List of the form [xmin, xmax, ymin, ymax, zmin, zmax]
    voxelSize : float
        Size of the voxels in Angstrom

    Returns
    -------
    grid : numpy.ndarray
        Grid of shape (nx, ny, nz)
    box_N : numpy.ndarray
        Number of voxels in each dimension

    """
    # increase grid by 0.5 to sample everything
    xrange = np.arange(boundingBox[0][0], boundingBox[1][0] + 0.5, step=voxelSize)
    yrange = np.arange(boundingBox[0][1], boundingBox[1][1] + 0.5, step=voxelSize)
    zrange = np.arange(boundingBox[0][2], boundingBox[1][2] + 0.5, step=voxelSize)

    gridpoints = np.zeros((xrange.shape[0] * yrange.shape[0] * zrange.shape[0], 3))
    i = 0
    for x in xrange:
        for y in yrange:
            for z in zrange:
                gridpoints[i][0] = x
                gridpoints[i][1] = y
                gridpoints[i][2] = z
                i += 1
    return gridpoints, (xrange.shape[0], yrange.shape[0], zrange.shape[0])


def get_bb(points):
    """Return bounding box from a set of points (N,3)

    Parameters
    ----------
    points : numpy.ndarray
        Set of points (N,3)

    Returns
    -------
    boundingBox : list
        List of the form [xmin, xmax, ymin, ymax, zmin, zmax]

    """
    minx = np.min(points[:, 0])
    maxx = np.max(points[:, 0])

    miny = np.min(points[:, 1])
    maxy = np.max(points[:, 1])

    minz = np.min(points[:, 2])
    maxz = np.max(points[:, 2])
    bb = [[minx, miny, minz], [maxx, maxy, maxz]]
    return bb


def get_all_protein_resids(pdb_file):
    """Return all protein residues from a pdb file

    Parameters
    ----------
    pdb_file : str
        Path to pdb file

    Returns
    -------
    resids : numpy.ndarray
        Array of protein resids old -> new

    """
    try:
        prot = Molecule(pdb_file)
    except:
        exit("could not read file")
    prot.filter("protein and not hydrogen")
    # mapping = prot.renumberResidues(returnMapping=True)
    return prot.get("index", "protein and name CA")


def get_all_metalbinding_resids(pdb_file):
    """Return all metal binding residues from a pdb file

    Parameters
    ----------
    pdb_file : str
        Path to pdb file

    Returns
    -------
    resids : numpy.ndarray
        id of resids that are metal binding

    """
    try:
        prot = Molecule(pdb_file)
    except:
        exit("could not read file")
    prot.filter("protein and not hydrogen")
    return prot.get(
        "index",
        sel="name CA and resname HIS HID HIE HIP CYS CYX GLU GLH GLN ASP ASH ASN GLN MET",
    )


def compute_average_p_fast(point, cutoff=0.25):
    """Using KDTree find the closest gridpoints

    Parameters
    ----------
    point : numpy.ndarray
        Point of shape (3,)
    cutoff : float
        Cutoff distance in Angstrom

    Returns
    -------
    average_p : numpy.ndarray
        Average probability of shape (1,)"""
    p = 0
    nearest_neighbors, indices = tree.query(
        point, k=20, distance_upper_bound=cutoff, workers=1
    )
    if np.min(nearest_neighbors) != np.inf:
        p = np.mean(output_v[indices[nearest_neighbors != np.inf]])
    return p


def get_probability_mean(grid, prot_centers, pvalues):
    """Compute the mean probability of all gridpoints from the globalgrid based on the individual boxes

    Parameters
    ----------
    grid : numpy.ndarray
        Grid of shape (nx, ny, nz)
    prot_centers : numpy.ndarray
        Protein centers of shape (N,3)
    pvalues : numpy.ndarray
        Probability values of shape (N,1)

    Returns
    -------
    mean_p : numpy.ndarray
        Mean probability over grid of shape (nx, ny, nz)
    """
    global output_v
    output_v = pvalues
    global prot_v
    prot_v = prot_centers
    cpuCount = multiprocessing.cpu_count()

    global tree
    tree = KDTree(prot_v)
    p = Pool(cpuCount)
    results = p.map(compute_average_p_fast, grid)
    return np.array(results)


def write_cubefile(bb, pvalues, box_N, outname="Metal3D_pmap.cube", gridres=1):
    """Write a cube file from a probability map
    The cube specification from gaussian is used, distance are converted to bohr

    Parameters
    ----------
    bb : list
        List of the form [xmin, xmax, ymin, ymax, zmin, zmax]
    pvalues : numpy.ndarray
        Probability values of shape (nx, ny, nz)
    box_N : tuple
        Number of voxels in each dimension
    outname : str
        Name of the output file
    gridres:float
        Resolution of the grid used for writing the voxels

    """

    with open(outname, "w") as cube:
        cube.write(" Metal3D Cube File\n")
        cube.write(" Outer Loop: X, Middle Loop y, inner Loop z\n")

        angstromToBohr = 1.89
        cube.write(
            f"    1   {bb[0][0]*angstromToBohr: .6f}  {bb[0][1]*angstromToBohr: .6f}   {bb[0][2]*angstromToBohr: .6f}\n"
        )
        cube.write(
            f"{str(box_N[0]).rjust(5)}    {1.890000*gridres:.9f}    0.000000    0.000000\n"
        )
        cube.write(
            f"{str(box_N[1]).rjust(5)}    0.000000    {1.890000*gridres:.9f}    0.000000\n"
        )
        cube.write(
            f"{str(box_N[2]).rjust(5)}    0.000000    0.000000    {1.890000*gridres:.9f}\n"
        )
        cube.write("    1    1.000000    0.000000    0.000000    0.000000\n")

        o = pvalues.reshape(box_N)
        for x in range(box_N[0]):
            for y in range(box_N[1]):
                for z in range(box_N[2]):
                    cube.write(f" {o[x][y][z]: .5E}")
                    if z % 6 == 5:
                        cube.write("\n")
                cube.write("\n")


def find_unique_sites(
    pvalues, grid, writeprobes=False, probefile="probes.pdb", threshold=7, p=0.1
):
    """The probability voxels are points and the voxel clouds may contain multiple metals
    This function finds the unique sites and returns the coordinates of the unique sites.
    It uses the AgglomerativeClustering algorithm to find the unique sites.
    The threshold is the maximum distance between two points in the same cluster it can be changed to get more metal points."""

    points = grid[pvalues > p]
    point_p = pvalues[pvalues > p]
    if len(points) == 0:
        print("no points available for clustering")
        return None
    clustering = AgglomerativeClustering(
        n_clusters=None, linkage="complete", distance_threshold=threshold
    ).fit(points)
    print(f"p={p}, n(metals):", clustering.n_clusters_)
    sites = []
    for i in range(clustering.n_clusters_):
        c_points = points[clustering.labels_ == i]
        c_points_p = point_p[clustering.labels_ == i]

        # compute center of probabilities as COM
        CM = np.average(c_points, axis=0, weights=c_points_p)
        position = CM
        sites.append((position, np.max(c_points_p)))
    if writeprobes:
        print(f"writing probes to {probefile}")
        with open(probefile, "w") as f:
            for i, site in enumerate(sites):
                f.write(
                    f"HETATM  {i+1:3} ZN    ZN A  1    {site[0][0]: 8.3f}{site[0][1]: 8.3f}{site[0][2]: 8.3f}  {site[1]:.2f}  0.0           ZN2+\n"
                )


def show_map(
    pdb,
    id,
    p=0.5,
    show_sticks_all=False,
    show_sticks_metalbinding=True,
    show_probes=True,
    show_pdb_metals=True,
):
    """Show the probability map of a protein and the protein using py3Dmol

    Parameters
    ----------
    pdb : str
        path of the pdb file
    id : str
        Name of the pdb file used as identifier for the probes and cube files
    p : float
        Isovalue of the probability map
    show_sticks_all : bool
        If True, show all the residues as sticks in the protein
    show_sticks_metalbinding : bool
        If True, show the residues that are metal binding as sticks in the protein

    Returns
    -------
    view : py3Dmol.view
        The view of the protein and the probability map
    """
    view = py3Dmol.view(width=1000, height=800)

    view.addModel(open(pdb, "r").read(), "pdb")
    if show_probes:
        view.addModel(open("probes_" + id + ".pdb", "r").read(), "pdb")
        probes = open("probes_" + id + ".pdb", "r").readlines()
        probabilites = [float(p[55:60]) for p in probes]  # read p from occupancy column
        colors = {}
        # use different colors for the probabilities
        for i, x in enumerate(probabilites):
            colors[i] = "#%02x%02x%02x" % (0, 0, int(x * 255 / 2))

    # check file size of cube file
    if os.path.getsize("metal_" + id + ".cube") < 10000000:
        view.addVolumetricData(
            open("metal_" + id + ".cube", "r").read(),
            "cube",
            {"isoval": p, "color": "blue", "opacity": 0.75},
        )
    else:
        print("cube file is too big, not showing the map")
        print(
            "Use smaller grid resolution or choose custom residues in close proximity to reduce file size. \n Alternatively, view full size cube file in VMD or UCSF ChimeraX."
        )

    view.zoomTo()
    view.setBackgroundColor("white")
    view.setStyle({}, {"cartoon": {"color": "gray"}})
    if show_sticks_all:
        view.setStyle({}, {"stick": {}, "cartoon": {"color": "gray"}})
    if show_pdb_metals:
        view.getModel(0).setStyle({"resn": "ZN"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "CA"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "CU"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "HG"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "MG"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "FE"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "MN"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "NI"}, {"sphere": {}})
        view.getModel(0).setStyle({"resn": "MB"}, {"sphere": {}})

    if show_probes:
        view.getModel(1).setStyle(
            {"resn": "ZN"},
            {"sphere": {"colorscheme": {"prop": "index", "map": colors}}},
        )

    # add hoverable labels for the residues and the predicted metals
    # two callbacks are needed, one for the residues and one for the metals
    # the metal one displays the probability
    view.getModel(0).setHoverable(
        {},
        True,
        """function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel(atom.resn+atom.resi+":"+atom.atom,{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                    }}""",
        """function(atom,viewer) { 
                    if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                    }""",
    )
    view.getModel(1).setHoverable(
        {},
        True,
        """function(atom,viewer,event,container) {
                    if(!atom.label) {
                        atom.label = viewer.addLabel('ZN p='+atom.pdbline.substring(55,60),{position: atom, backgroundColor: 'mintcream', fontColor:'black'});
                    }}""",
        """function(atom,viewer) { 
                    if(atom.label) {
                        viewer.removeLabel(atom.label);
                        delete atom.label;
                    }
                    }""",
    )
    if show_sticks_metalbinding:
        view.setStyle({"resn": "HIS"}, {"stick": {}, "cartoon": {"color": "gray"}})
        view.setStyle({"resn": "ASP"}, {"stick": {}, "cartoon": {"color": "gray"}})
        view.setStyle({"resn": "GLU"}, {"stick": {}, "cartoon": {"color": "gray"}})
        view.setStyle({"resn": "CYS"}, {"stick": {}, "cartoon": {"color": "gray"}})

    return view.show()


def show_probes(probefile):
    """visualizes the probefile in vmd

    Parameters
    ----------
    probefile: str
        path to probefile

    """
    try:
        probes = Molecule(probefile)
    except:
        exit("could not read probefile")
    probes.view(style="VDW", color="Occupancy", viewer="vmd")


def maxprobability(pvalues, grid, pdb, label):
    """Returns the point in the grid with the highest probability"""

    print(
        f"Max p={np.max(pvalues)} xyz: {grid[np.argmax(pvalues)][0]: .4f} {grid[np.argmax(pvalues)][1]: .4f} {grid[np.argmax(pvalues)][2]: .4f}"
    )
    if os.path.isfile("maxp_" + os.path.basename(label) + ".csv"):
        with open("maxp_" + os.path.basename(label) + ".csv", "a") as fp:
            fp.write(
                f"{pdb},{np.max(pvalues):.4f}, {grid[np.argmax(pvalues)][0]: .4f}, {grid[np.argmax(pvalues)][1]: .4f}, {grid[np.argmax(pvalues)][2]: .4f}\n"
            )
    else:
        with open("maxp_" + os.path.basename(label) + ".csv", "w") as fp:
            fp.write("pdb,p,x,y,z\n")
            fp.write(
                f"{pdb},{np.max(pvalues):.4f}, {grid[np.argmax(pvalues)][0]: .4f}, {grid[np.argmax(pvalues)][1]: .4f}, {grid[np.argmax(pvalues)][2]: .4f}\n"
            )
