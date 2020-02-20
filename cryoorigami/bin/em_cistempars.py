#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Date    : 2018-11-15 11:00:21
# @Author  : Your Name (you@example.org)
# @Link    : http://example.org
# @Version : $Id$

import os
import sys
import argparse
import cryoorigami.utilities as util


def main():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument("-fsc50",    "--fsc50",       type=float, help="FSC50 from fsc curve. Get closest value above 0.5", default=None)
    parser.add_argument("-plimit",   "--plimit",      type=float, help="Previous resolution limit", default=20.0)
    parser.add_argument("-ppercent", "--ppercent",    type=float, help="Previous percentage of particles", default=0)
    parser.add_argument("-n",        "--nparticles",  type=int,   help="Total number of particles")
    parser.add_argument("-r",        "--radius",      type=float, help="Radius of the spherical mask applied to 3D map")
    parser.add_argument("-k",        "--numclasses",  type=int,   help="Number of classes", default=1)

    args = parser.parse_args()

    # Prepare args dict
    args_dict = {'fsc50':       args.fsc50,
                 'plimit':      args.plimit,
                 'ppercent':    args.ppercent,
                 'nparticles':  args.nparticles,
                 'radius':      args.radius,
                 'numclasses':  args.numclasses
                 }

    new_limit, p_betterRes, p_worseRes = util.estimate_cistem_params(args_dict['fsc50'], 
                                                                     args_dict['ppercent'],
                                                                     args_dict['plimit'],
                                                                     args_dict['nparticles'],
                                                                     args_dict['radius'],
                                                                     args_dict['numclasses'])
    # Print parameters
    if new_limit is not None:
        print('New resolution limit: %.2f' % (new_limit))
    print('New p-particles for better-res: %.2f' % (p_betterRes))
    print('New p-particles for worse-res: %.2f' % (p_worseRes))


if __name__ == "__main__":
    main()
