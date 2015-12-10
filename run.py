#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import re
import itertools
from collections import OrderedDict

import execu
import torque
#########1#########2#########3#########4#########5#########6#########7#########


def args_latest():
    const = ['-N16384']
    params = OrderedDict()
    params.update(D=[2, 3])
    params.update(C=['neumann', 'moore', 'hex'])
    params.update(P=['push', 'pushfill', 'walkfill'])
    params.update(S=['random', 'poisson', 'even'])
    return [const + x + [make_label(x)] for x in product(params)]


def sequential(params):
    ret = []
    for (key, vals) in params.items():
        var_args = []
        for value in vals:
            if len(key) > 1:
                var_args.append('--{}={}'.format(key, value))
            else:
                var_args.append('-{}{}'.format(key, value))
            ret.append(var_args)
    return ret


def product(params):
    ret = []
    for vals in itertools.product(*params.values()):
        var_args = []
        for (key, value) in zip(params.keys(), vals):
            if len(key) > 1:
                var_args.append('--{}={}'.format(key, value))
            else:
                var_args.append('-{}{}'.format(key, value))
        ret.append(var_args)
    return ret


def make_label(var_args):
    label = '_'.join([s.lstrip('-') for s in var_args])
    return '--label=' + re.sub('[^\w\._]+', '_', label)


#########1#########2#########3#########4#########5#########6#########7#########
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-v', '--verbose', action='store_true')
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-Q', '--torque', action='store_true')
    parser.add_argument('-C', '--directory', default=os.getcwd())
    parser.add_argument('-q', '--queue',
                        choices=['low', 'batch', 'high'], default='batch')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    parser.add_argument('outfile', nargs='?', default=sys.stdout)
    args = parser.parse_args()

    program = os.path.join(os.path.dirname(__file__), 'a.out')
    constargs = [program]
    constargs.append('--top_dir=' + args.directory)

    args_list = args_latest()
    commands = [constargs + x for x in args_list] * args.repeat

    if args.torque:
        qargs = dict()
        qargs['-q'] = args.queue
        rex = re.compile('--label=(\S+)')
        for cmd in commands:
            qargs['-N'] = rex.search(' '.join(cmd)).group(1)
            torque.qsub(cmd, args.dry_run, ppn=1, **qargs)
    else:
        pool = execu.Pool()
        for cmd in commands:
            strcmd = ' '.join(cmd)
            print(strcmd)
            if not args.dry_run:
                pool.apply_async(cmd, name=' '.join(cmd))

    print('End of ' + __file__)
