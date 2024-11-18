from configparser import ConfigParser
from difflib import ndiff
import hashlib
from datadispatch.access import select

# TODO: 
SNAKEFOOD = '/home/erboone/pipeline/snakefood.ini'
DEFAULT_SNAKEFOOD = '/home/erboone/pipeline/snakefood.default.ini'

def order_snakefood(
        section_key:str=None, 
        alt_name:str=None,
        default:bool=False) -> dict | ConfigParser:
    
    config = ConfigParser()
    if default:
        config.read(DEFAULT_SNAKEFOOD)
    elif alt_name:
        config.read(alt_name)
    else:
        config.read(SNAKEFOOD)

    if section_key:
        return config[section_key]
    else:
        return config

def hash_strli(l:list[str], num_chars:int=8) -> str:
    b_l = bytes('_'.join(l), encoding='utf-8')
    h = hashlib.sha256(b_l).hexdigest()[-8:]

    return h


# Hash function
def hash_parameters(ini_section_key:str, alt_name:str=None) -> str:
    parameters = (
        dict(sorted(order_snakefood(ini_section_key, alt_name).items())),
        dict(sorted(order_snakefood(ini_section_key, alt_name, default=True).items()))
    )

    active_par, default_par = parameters
    if active_par.keys() != default_par.keys():
        raise RuntimeError('Mismatch in keys bewteen default and active snakefood files\nActive: {}\nDefault: {}'.
                           format(*[active_par.keys(), default_par.keys()]))
    
    diff = ndiff(
        list(default_par.values()),
        list(active_par.values())
    )
    
    changes = []
    for i in diff:
        change = i[0]
        text = i[2:]
        if change == '-':
            continue
        elif change == '+':
            changes.append(text)
        elif change == ' ':
            changes.append('')
    
    return hash_strli(changes)
    

# def translate 

if __name__ == "__main__":
    from A_CellSeg import *

    A1_Cellpose_out = 'M182_.checkpoint'
    Cellpose('202411080959_20241108M182MusLivSenTS2_VMSC02201', 'M182_0fb68fde.checkpoint')