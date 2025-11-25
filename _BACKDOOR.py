# from scripts import A_CellSeg

# A_CellSeg.Cellpose('202403201153_20240320M151MUSBrSenTS2_VMSC02201', 'test.out', (100, 150))

from datetime import datetime
from subprocess import run
import re
from pathlib import Path
import os
from glob import glob

# from scripts import _1_
# import scripts.MERSCOPE._1_Segmentation as _1_ 
import scripts.MERSCOPE._2_Preproc as _2_
import scripts.MERSCOPE._3_Annotation as _3_

from scripts.meta import *
from datadispatch.access import select
from mftools.mfexperiment import MerscopeExperiment, SmallMerscopeExperiment
from string import Template

hashes: dict[str, object] = {}
all_config = order_snakefood()
check_RunInfo(all_config)
print('all_config', all_config.keys())
for stepname in all_config.keys():
    h = hash_parameters(stepname)  # in scripts.meta
    hashes[stepname] = h

commit = get_commit_hash("TODO: this parameter does nothing")


def _1_Merquaco(name,
        input,
        outfile,
        hashes,
        commit,):
    pass

def cardioids(code):
    BASE = '/data/erboone/cardioids'
    CODE = code
    _3_.A4_HarmAnnotation([f'{BASE}/{CODE}.h5ad'], f'{BASE}/{CODE}.ann.h5ad', None, None)
    _3_.QC_2_postanno(f'{BASE}/{CODE}.ann.h5ad', f"{BASE}/_output/{CODE}/stats.csv", None, None)

def proseg():
    BASE = '/data/erboone/proseg/'
    NAMES = ['202405211218_HumanBrain-FF-UCSD-CP1403-PU-RC_VMSC03501_new', 
             '202409091924_240904JHCX220206G8-HuBrain-CX22-37C-LH_VMSC12502',
             '202409111637_240904JHCX22021L4G--HuBrain-CX22-LH_VMSC11802',
             '202411231816_241118JHCX2303EK3-BICAN-VS262-2006-LH_VMSC18710']
    for name in NAMES:
        _3_.A4_HarmAnnotation([f'{BASE}{name}/adata.h5ad'], f'{BASE}{name}/adata.anno.h5ad', None, None)
        _3_.QC_2_postanno(f'{BASE}{name}/adata.anno.h5ad', f"{BASE}{name}/adata.anno/stats.csv", None, None)

def default():
    BASE = '/home/erboone/pipeline/_output/B3_Annotation'
    NAMES = ['/home/erboone/pipeline/_output/B3_Annotation/Ren29.efb3ba87.51c55064.h5ad']
    for name in [Path(n) for n in NAMES]:
        # _3_.A4_HarmAnnotation([f'{BASE}{name}/adata.h5ad'], f'{BASE}{name}/adata.anno.h5ad', None, None)
        # _3_.QC_2_postanno(f'{BASE}/{name}.h5ad', f"{BASE}/{name}/stats.csv", None, None)
        _3_.QC_2_postanno(name, name.parent /name.stem / "stats.csv", None, None)


def sennet():
    Start = '/home/erboone/SenNet/jeolness_cp_20250516/250512_SenTMLiv_126to194_pre.h5ad'
    out = '/home/erboone/SenNet/annotations'
    # DANGER this is not correct
    _2_.B1_SaveRawScanpy([Start], f'{out}/126to194.anno.h5ad', None, None, None)
    _2_.B2_Preprocessing([Start], f'{out}/126to194.anno.h5ad', None, None, None)
    _3_.A4_HarmAnnotation([Start], f'{out}/126to194.anno.h5ad', None, None)
    _3_.QC_2_postanno(f'{out}/adata.anno.h5ad', f"{out}/adata.anno/stats.csv", None, None)



if __name__ == '__main__':
    # res = select('Experiment', where="Experiment.name == '202502281312_20250228M189BICANRen24_VMSC32110'")

    # res = res.pop()
    # _1_.A1_Cellpose(res.name, 0, 'thisWorks.txt', 0, 0, 0)

    # _A1_._2_SaveRawScanpy(res.name, 0, output='test', hashes={'A11_Cellpose':'test_hash'},commit='test_commit')
    # _A2_._1_HarmAnnotation(["/home/erboone/pipeline/_output/A1_Segmentation/T001_a436536e.h5ad"], 'test)
    # a = MerscopeExperiment(res.rootdir, res.name).create_scanpy_object()
    # print(a.X.dtype)

    # _2_.B1_SaveRawScanpy('202505311513_20250531M200BICANRen32_VMSC32010', None, '/home/erboone/pipeline/_output/B2_Preproc/Ren32.h5ad', None, '5')
    #_2_.B2_Preprocessing(#'/home/erboone/pipeline/_output/A1_Segmentation/T001_0558562c.h5ad'
    #     res.name, '/home/erboone/pipeline/_output/B2_Preproc/T001.h5ad', '/home/erboone/pipeline/_output/B2_Preproc/T001.filt.h5ad', None, None)
    
    # _2_.QC_1_postsegqc('/home/erboone/pipeline/_output/B2_Preproc/Ren24_da4470d6.filt.h5ad', '/home/erboone/pipeline/_output/B2_Preproc/Ren24_da4470d6.QC.png', None, None)
    # ultra_mfadata = SmallMerscopeExperiment("/mnt/merfish15/MERSCOPE/vizgen_MERSCOPE/", "202411231816_241118JHCX2303EK3-BICAN-VS262-2006-LH_VMSC18710/region_R1/").create_scanpy_object()
    # ultra_mfadata.obs['region'] = "VIZ008R1"
    # ultra_mfadata.obs_names = ultra_mfadata.obs_names.astype(str)
    # ultra_mfadata.layers['counts'] = ultra_mfadata.X
    # ultra_mfadata.write('tempVIZ008R1.h5ad')



    # _3_.A4_HarmAnnotation(['/data/erboone/proseg/Ren26/adata.h5ad'], f'/data/erboone/proseg/Ren26/adata.anno.h5ad', None, None)
    # _3_.QC_2_postanno(f'/data/erboone/proseg/Ren26/adata.anno.h5ad', f"/data/erboone/proseg/Ren26/adata.anno/stats.csv", None, None)

    # _3_.A4_HarmAnnotation(['/data/erboone/proseg/Ren27/adata.h5ad'], f'/data/erboone/proseg/Ren27/adata.anno.h5ad', None, None)
    # _3_.QC_2_postanno(f'/data/erboone/proseg/Ren27/adata.anno.h5ad', f"/data/erboone/proseg/Ren27/adata.anno/stats.csv", None, None)
    # proseg()

    # _3_.QC_2_postanno(f'/home/erboone/pipeline/_output/B3_Annotation/Ren31.efb3ba87.16264781.h5ad', f"/home/erboone/pipeline/_output/B3_Annotation/Ren31.efb3ba87.16264781/stats.csv", None, None)
    # Viz1_putamen = SmallMerscopeExperiment("/mnt/merfish15/MERSCOPE/vizgen_MERSCOPE/", "202405211218_HumanBrain-FF-UCSD-CP1403-PU-RC_VMSC03501_new/region_R1/")
    # Viz4_putamen = SmallMerscopeExperiment("/mnt/merfish15/MERSCOPE/vizgen_MERSCOPE/", "202409091924_240904JHCX220206G8-HuBrain-CX22-37C-LH_VMSC12502/region_R1/")
    # Viz5_putamen = SmallMerscopeExperiment("/mnt/merfish15/MERSCOPE/vizgen_MERSCOPE/", "202409111637_240904JHCX22021L4G--HuBrain-CX22-LH_VMSC11802/region_R1/")
    # Viz8_putamen = SmallMerscopeExperiment("/mnt/merfish15/MERSCOPE/vizgen_MERSCOPE/", "202411231816_241118JHCX2303EK3-BICAN-VS262-2006-LH_VMSC18710/region_R2/")

    # experiments = [
    #     ("viz1": Viz1_putamen),
    #     ("viz4": Viz4_putamen),
    #     ("viz5": Viz5_putamen),
    #     ("viz8": Viz8_putamen)
    # ]
    # _3_.A4_HarmAnnotation([f'/home/erboone/pipeline/_output/B2_Preproc/202405211218_HumanBrain-FF-UCSD-CP1403-PU-RC_VMSC03501_new.efb3ba87.filt.h5ad'], f'/home/erboone/pipeline/_output/B3_Annotation/202405211218_HumanBrain-FF-UCSD-CP1403-PU-RC_VMSC03501_new.efb3ba87.37e5737d.h5ad', None, None)
    # _3_.QC_2_postanno(f'{BASE}{name}/adata.anno.h5ad', f"{BASE}{name}/adata.anno/stats.csv", None, None)
    # _3_.QC_2_postanno(f'/home/erboone/pipeline/_output/B3_Annotation/Ren29.efb3ba87.16264781.h5ad', f"/home/erboone/pipeline/_output/B3_Annotation/Ren29.efb3ba87.16264781/stats.csv", None, None)

    # default()
    # _3_.QC_2_postanno(f'/home/erboone/pipeline/_output/B3_Annotation/Ren29.efb3ba87.51c55064.h5ad', '/home/erboone/pipeline/_output/B3_Annotation/Ren29.efb3ba87.51c55064/stats.csv', None, None)
    # paths = [
    #     ("/home/erboone/pipeline/_output/B3_Annotation/Ren36.efb3ba87.51c55064.h5ad", ""),
    # ]
    # _3_.A4_HarmAnnotation([f'/home/erboone/BICAN/temp/Ren14.h5ad'], f'/home/erboone/BICAN/temp/Ren14.anno.h5ad', None, None)
    # _3_.QC_2_postanno(f'/home/erboone/BICAN/temp/Ren14.anno.h5ad', f'/home/erboone/BICAN/temp/Ren14.anno/stats', None, None)
    # from scripts.meta import read_ref_table
    # read_ref_table('Pu')
    # _2_.C2_Preprocessing("/home/erboone/pipeline/_output/C2_Preproc/C2.checkpoint", None, None, None)
    # _2_.QC_1_postsegqc()
    cardioids()

    

