import argparse
from sscpol.fitter import SSC_Fitter

parser = argparse.ArgumentParser()
parser.add_argument('--method', type=str, choices=["standard","direct", "ps"], default="standard",
                    help='What optimizer to use')
parser.add_argument('--blazar', type=str, choices=["J2011","TXS", "S5low", "S5flare"], default="J2011",
                    help='Which blazar to fit')
parser.add_argument('--seed', type=int, default=42,
                    help='Random seed for magnetic field zones')
parser.add_argument('--nblocks', type=int, choices=[1,7,19,37,61,91,127], default=1,
                    help='Number of magnetic field zones')
args = parser.parse_args()

def main():
    print("BLAZAR:  ", args.blazar)
    fitter = SSC_Fitter(args.blazar, method=args.method, seed=args.seed, nblocks=args.nblocks)
    res = fitter.run()
    print("FIT RESULT: ", res)

if __name__ == "__main__":
    main()