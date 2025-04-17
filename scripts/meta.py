import os
from configparser import ConfigParser
from datadispatch.access import update, add
from datadispatch.orm import ParamLog

# Constants
SNAKEFOOD = '/home/erboone/pipeline/snakefood.ini'
DEFAULT_SNAKEFOOD = '/home/erboone/pipeline/snakefood.default.ini'
_output = f'{os.getenv("HOME")}/pipeline/_output'
OUTPUT = _output
RUNINFO_KEY = 'Run Info'
REQ_RUNINFO_VALUES = ['description', 'where_description', 'where']

# File formats (to reduce typos)
S_CHECKPOINT = 'checkpoint'
S_H5AD = 'h5ad'
S_IMAGE = 'png'
S_CSV = 'csv'

# Error strings
ERR_NO_RUNINFO_SECTION = "Every snakefood file requires a section 'Run Info'"
ERR_NO_RUNINFO_VALUE = ""
ERR_MISSING_RUNINFO_KEY = \
    f"Every snakefood file requires a section 'Run Info' with keys {REQ_RUNINFO_VALUES} and corresponding string values. \nFound:\n{{}}"
ERR_MYSTERY_CONFIG_INPUT = "Make sure all 'Run Info' values are filled, with strings preferably."

def order_snakefood(
        section_key:str=None, 
        alt_name:str=None,
        default:bool=False) -> dict:
    
    config = ConfigParser()
    if default:
        config.read(DEFAULT_SNAKEFOOD)
    elif alt_name:
        config.read(alt_name)
    else:
        config.read(SNAKEFOOD)

    if not config.sections():
        raise RuntimeError(f"config file load failed, check '{DEFAULT_SNAKEFOOD}', '{alt_name}', or '{SNAKEFOOD}'")
    

    if section_key:
        return config[section_key]
    else:
        return config

def hash_strli(l:list[str], num_chars:int=8) -> str:
    b_l = bytes('_'.join(l), encoding='utf-8')
    h = hashlib.sha256(b_l).hexdigest()[-8:]

    return h


# Hash function
from difflib import ndiff
import hashlib
def _old_hash_parameters(ini_section_key:str, alt_name:str=None) -> str:
    # Depreciated
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


def hash_parameters(ini_section_key:str, alt_name:str=None) -> str:
    parameters = (
        dict(sorted(order_snakefood(ini_section_key, alt_name).items())),
        dict(sorted(order_snakefood(ini_section_key, alt_name, default=True).items()))
    )
    active_par, default_par = parameters

    if active_par.keys() != default_par.keys():
        raise RuntimeError('Mismatch in keys bewteen default and active snakefood files\nActive: {}\nDefault: {}'.
                           format(*[active_par.keys(), default_par.keys()]))
    
    param_strl = [str(p) for p in active_par.values()]

    return hash_strli(param_strl)


import subprocess as sub
def get_commit_hash(git_folder_path):
    _pipe = sub.Popen(['git', 'log', '-1'], stdout=sub.PIPE)
    _pipe = sub.Popen(['head', '-1'], stdout=sub.PIPE, stdin=_pipe.stdout)
    out_bytes = sub.check_output(['awk',  '{print $2}'], stdin=_pipe.stdout)
    commit = out_bytes.decode('utf-8').strip()
    return commit

def branch_up_to_date():
    _temp = sub.Popen(['git', 'status', '--short'], stdout=sub.PIPE)
    out_bytes = sub.check_output(['awk',  '{print $1}'], stdin=_temp.stdout)
    return 'M' in set(out_bytes.decode('utf-8').strip().split('\n'))



from string import Template
from itertools import product
# for easy alteration; base_T = base template can be changed here or in func
# call
base_T = '{_o}/{_step}/{_file}'
def assemble_target(
        template:dict,
        hashes:dict,
        format:str,
        wildcards:tuple[str, list],
        _base:str='{_o}/{_step}/{_file}'
    ):

    # NOTE: Excecution pattern is 
    ### ex.) step description -> Variable example ###
    # 1.) Define for step -> A1_T:Template
    # 2.) Fill hashes
    # 3.) Fill formats 
    # 4.) Full wildcards
    
    # validate params
    req_keys = ['_o', '_step', '_file']
    if any([k not in req_keys for k in template.keys()]):
        raise RuntimeError(f"template parameter requires {req_keys} as keys")

    if '$' not in template['_file']:
        raise RuntimeError(f"template['_file'] parameter requires a $-based substitution (i.e. ${{identifier}}), see python documentation: https://docs.python.org/3.4/library/string.html#template-strings")

    # step 1: -----------------------------
    T = Template(base_T.format(**template))

    # step 2: -----------------------------
    H = Template(T.safe_substitute(hashes))

    # step 3: -----------------------------
    F = Template(H.safe_substitute({'form':format}))
    F_return = F.safe_substitute({
        ident:'{' + f'{ident}' + '}'
        for ident in F.get_identifiers()
    })

    # step 4: -----------------------------

    # keep this in case product of multiple WC is required
    # param: wildcards:dict[list],
    # test = [zip([k] * len(v), v) for k, v in wildcards.items()]
    if not wildcards:
        W_return = F_return
    else:
        W_return = [
            F.safe_substitute({
                wildcards[0]:i
            }) for i in wildcards[1]]
    

    return W_return , F_return


def log_parameters(hashes, config):
    new_entrys = []
    for k, v in hashes.items():
        if k in ['SNAKEMAKE_CALL', ]:
            pass
        
        if k.isupper():
            new_entry = ParamLog(
                runname=hashes['SNAKEMAKE_CALL'],
                step=k,
                hash='--------',
                config=str(v)
            )
        else:
            new_entry = ParamLog(
                runname=hashes['SNAKEMAKE_CALL'],
                step=k,
                hash=v,
                config=str(dict(config[k]))
            )
        
        new_entrys.append(new_entry)

    add(new_entrys)

def mark_success(dt:str):
    update(ParamLog, 
           {'success':True}, 
           f"ParamLog.runname == '{dt}'")

def check_RunInfo(config):
    
    try:
        ri = config[RUNINFO_KEY]
    except KeyError as e:
        raise KeyError(ERR_NO_RUNINFO_SECTION) from e

    for key in REQ_RUNINFO_VALUES:
        try:
            if ri[key] == '':
                raise RuntimeError
            
        except (KeyError, RuntimeError) as e:
            kvs = dict(ri).items()
            s = ''.join([f'\t{k}:\t{v}\n' for k,v in kvs])
            raise RuntimeError(ERR_MISSING_RUNINFO_KEY.format(s)) from e
        
        except NotImplementedError as e:
            raise RuntimeError(ERR_MYSTERY_CONFIG_INPUT) from e


if __name__ == "__main__":
    print(assemble_target(
        template={
            '_o': _output,
            '_step': 'A1_Segmentation',
            '_file':'${exp_id}_${fake_wc}_${a11}.${form}',
        },
        hashes={
            'a11':'#narhsh'
        },
        format='checkpoint',
        wildcards=('exp_id',['1', '2', '3'])
    ))
