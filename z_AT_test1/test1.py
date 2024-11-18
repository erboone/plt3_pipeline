#!/usr/bin/env python

from ..scripts.helper import *

from configparser import ConfigParser

config = ConfigParser()
config.read('../snakefood.ini')
CONFIG = config['test1']

ONE = CONFIG['1']
TWO = CONFIG['2']
THREE = CONFIG['3']
FOUR = CONFIG['4']


hash = hash_parameters((ONE, TWO, THREE, FOUR))
filename = '_output/{}'.format(*[ONE, TWO, THREE, FOUR])


with open(filename, 'w') as out:
    out.write('{}'.format())
