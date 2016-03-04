#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
"""
import sys
import os
import glob
import re
import itertools
import datetime
import subprocess
import multiprocessing
from collections import OrderedDict

import execu
import torque
#########1#########2#########3#########4#########5#########6#########7#########


def args_latest():
    return args_k()


def args_all():
    const = ['-N16384']
    params = OrderedDict()
    params.update(D=[2, 3])
    params.update(C=['neumann', 'moore', 'hex'])
    params.update(P=['push', 'pushn', 'pushne', 'fillpush', 'fill', 'empty'])
    params.update(S=['random', 'even'])
    return [const + x + ['--out_dir=' + make_outdir(x,i)] for i,x in enumerate(product(params))]


def args_k():
    const = ['-v']
    params = OrderedDict()
    params['k'] = [10 ** x for x in range(6)]
    print(params)
#    return sequential(params)
    return [const + [x, '--out_dir=' + make_outdir([x])] for x in sequential(params)]


def sequential(params):
    for (key, vals) in params.items():
        for value in vals:
            if len(key) > 1:
                yield '--{}={}'.format(key, value)
            else:
                yield '-{}{}'.format(key, value)


def product(params):
    for vals in itertools.product(*params.values()):
        var_args = []
        for (key, value) in zip(params.keys(), vals):
            if len(key) > 1:
                var_args.append('--{}={}'.format(key, value))
            else:
                var_args.append('-{}{}'.format(key, value))
        yield var_args


def make_outdir(var_args=[], i=0):
    prefix = 'tumopp'
    label = '_'.join([s.lstrip('-') for s in var_args])
    now = datetime.datetime.now().strftime('%Y%m%d-%H%M')
    pid = '{}-{:04}'.format(os.getpid(), i)
    return '_'.join([prefix, label, now, pid])


#########1#########2#########3#########4#########5#########6#########7#########
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-j', '--threads', type=int,
                        default=multiprocessing.cpu_count())
    parser.add_argument('-B', '--batch', action='store_const', const='batch')
    parser.add_argument('-Q', '--torque', action='store_const', const='torque', dest='batch')
    parser.add_argument('-q', '--queue',
                        choices=['low', 'batch', 'high'], default='batch')
    parser.add_argument('-r', '--repeat', type=int, default=1)
    (args, rest) = parser.parse_known_args()

    project = os.path.dirname(__file__)
    program = os.path.join(project, 'a.out')
    postproc = os.path.join(project, 'post.R')
    constargs = [program]

    if not args.batch:
        outdir = make_outdir(rest)
        cmd = constargs + rest + ['-v', '--out_dir=' + outdir]
        print(' '.join(cmd))
        if not args.dry_run:
            subprocess.call(cmd)
            for od in glob.glob(outdir + '*'):
                subprocess.call([postproc, od])
        exit()

    args_list = args_latest()
    commands = [constargs + x for x in args_list] * args.repeat

    rex = re.compile(r'--out_dir=(\S+)')
    if args.batch == 'torque':
        qargs = dict()
        qargs['-q'] = args.queue
        for cmd in commands:
            strcmd = ' '.join(cmd)
            outdir = rex.search(strcmd).group(1)
            qargs['-N'] = outdir
            torque.qsub(cmd, args.dry_run, ppn=1, **qargs)
    elif args.batch == 'batch':
        pool = execu.Pool(args.threads)
        for cmd in commands:
            strcmd = ' '.join(cmd)
            outdir = rex.search(strcmd).group(1)
            print(strcmd)
            if not args.dry_run:
                pool.apply_async(cmd, name=outdir)

    print('End of ' + __file__)
