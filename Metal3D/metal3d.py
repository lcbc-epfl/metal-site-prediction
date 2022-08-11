#!/usr/bin/env python3

import glob
import argparse
import time

import torch
import torch.nn as nn
from utils.helpers import *

from utils.voxelization import processStructures as processStructures
from utils.model import Model
from moleculekit.tools.voxeldescriptors import getVoxelDescriptors, viewVoxelFeatures

from moleculekit.util import boundingBox

from scipy.spatial import distance

import warnings


if __name__ == "__main__":
    start_time = time.time()

    parser = argparse.ArgumentParser(description="Inference using Metal3D")
    parser.add_argument(
        "--id",
        nargs="+",
        help="indexes of CA atoms to predict",
        type=int,
        default=[],
    )
    parser.add_argument(
        "--metalbinding",
        help="uses all residues that have sidechains that can coordinate metals",
        action="store_true",
    )
    parser.add_argument("--writecube", help="write a cube file", action="store_true")
    parser.add_argument(
        "--cubefile", help="name of cube file", default="Metal3D_pred.cube"
    )
    parser.add_argument(
        "--no-clean", help="do not clean the pdb file", action="store_false"
    )
    parser.add_argument("--pdb", help="pdb, pdb.gz, cif or cif.gz file", required=True)
    parser.add_argument("--batch-size", help="batchsize", default=128, type=int)

    parser.add_argument("--writeprobes", help="write probe files", action="store_true")
    parser.add_argument(
        "--probefile",
        help="name of files with predicted metal positions ",
        default="probes_Metal3D.pdb",
    )
    parser.add_argument(
        "--probeprobabilities",
        nargs="+",
        help="probabilities to predict probes",
        type=float,
        default=[0.1, 0.2, 0.3, 0.5, 0.75],
    )
    parser.add_argument("--threshold", help="cluster threshold", type=float, default=7)
    parser.add_argument("--pthreshold", help="p threshold", type=float, default=0.10)
    parser.add_argument(
        "--maxp", help="print max probability point", action="store_true"
    )
    parser.add_argument("--label", help="label for the maxp file", default="metal3d")
    parser.add_argument(
        "--softexit",
        help="dont ask for conformation before exiting",
        action="store_true",
    )
    args = parser.parse_args()

    if len(args.id) == 0 and args.metalbinding == False:
        print("No resid passed, using whole protein")
        ids = get_all_protein_resids(args.pdb)
    elif len(args.id) == 0 and args.metalbinding:
        print("Using all residues that can bind metals with their sidechain")
        ids = get_all_metalbinding_resids(args.pdb)
    else:
        print("Using the following indexes", args.id)
        ids = args.id

    voxels, prot_centers, prot_N, prots = processStructures(args.pdb, ids)

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    model = Model()
    model.to(device)
    model.load_state_dict(torch.load(f"weights/metal_0.5A_v3_d0.2_16Abox.pth"))
    voxels.to(device)
    model.eval()
    outputs = torch.zeros([voxels.size()[0], 1, 32, 32, 32])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore")
        print(voxels.shape)
        for i in range(0, voxels.size()[0], args.batch_size):
            o = model(voxels[i : i + args.batch_size])
            outputs[i : i + args.batch_size] = o.cpu().detach()

    # process all predicted probabilities
    prot_v = np.vstack(prot_centers)
    output_v = outputs.flatten().numpy()

    bb = get_bb(prot_v)

    grid, box_N = create_grid_fromBB(bb)

    probability_values = get_probability_mean(grid, prot_v, output_v)

    if args.writecube:
        write_cubefile(
            bb,
            probability_values,
            box_N,
            outname=args.cubefile,
            gridres=1,
        )

    if not args.softexit:
        viewVoxelFeatures(probability_values, grid, box_N, featurenames=[f"p(metal)"])

    if not args.softexit:
        prots[0].view(style="NewCartoon", hold=True, viewer="vmd")
        prots[0].view(
            sel="resname HIS HID HIE HIP ASP ASH GLU GLH CYS CYX CYN and not hydrogen",
            style="Licorice",
            hold=False,
            viewer="vmd",
        )
    if args.writeprobes:
        find_unique_sites(
            probability_values,
            grid,
            writeprobes=args.writeprobes,
            probefile=args.probefile,
            threshold=args.threshold,
            p=args.pthreshold,
        )
        if not args.softexit:
            show_probes(probefile=args.probefile)

    if args.maxp:
        maxprobability(probability_values, grid, args.pdb, args.label)

    print("--- %s seconds ---" % (time.time() - start_time))

    if not args.softexit:
        input("Press Enter to end program...")
