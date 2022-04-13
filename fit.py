import argparse
from sscpol.fitter import SSC_Fitter

parser = argparse.ArgumentParser()
parser.add_argument('--method', type=str, choices=["standard","direct", "ps","cross_entropy"], default="standard",
                    help='What optimizer to use')
parser.add_argument('--blazar', type=str, choices=["J2011","TXS", "S5low", "S5flare", "J2011flare"], default="J2011",
                    help='Which blazar to fit')
parser.add_argument('--seed', type=int, default=500,
                    help='Number of random seeds for magnetic field zones')
parser.add_argument('--nblocks', type=int, choices=[1,7,19,37,61,91,127], default=1,
                    help='Number of magnetic field zones')
parser.add_argument('--nprocs', type=int, default=20, choices=[20,30],
                    help='Number of cpu processes')
parser.add_argument('--rand_gamma', type=int, default=0, choices=[0,1],
                    help='Whether individual zones have random maximum electron energies')
args = parser.parse_args()

def main():
    print("BLAZAR:  ", args.blazar)
    print("METHOD:  ", args.method)
    print("NBLOCKS:    ", args.nblocks)
    print("RAND_GAMMA: ", args.rand_gamma)
    fitter = SSC_Fitter(args.blazar, method=args.method, seed=args.seed, 
                        nblocks=args.nblocks, nprocs=args.nprocs, rand_gamma=args.rand_gamma)
    res = fitter.run()
    print("FIT RESULT: ", res)

if __name__ == "__main__":
    main()