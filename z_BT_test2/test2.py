#!/usr/bin/env python

from configparser import ConfigParser

config = ConfigParser()
config.read('../snakefood.ini')
CONFIG = config['test2']

ONE = CONFIG['1']
TWO = CONFIG['2']
THREE = CONFIG['3']
FOUR = CONFIG['4']

filename = '_output/{}{}{}{}'.format(*[ONE, TWO, THREE, FOUR])
with open(filename, 'w') as out:
    out.write('{}{}{}{}'.format(*['test1'.upper(), ONE, TWO, THREE, FOUR]))
