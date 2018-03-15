"""
FireWorks implementation for the computation of electronic structure for
molecules in the Metatlas database.

-----------    ----------    -----------
| Create  |    | Run    |    | Process |
| Orca    | => | Orca   | => | Output  |
| Input   |    | Calc   |    | File    |
-----------    ----------    -----------
"""
from fireworks import Firework
from metatlas import OBUFFOptimize, OrcaOptimize, \
    create_launchpad, make_df_with_smiles_only_from_csv


if __name__ == "__main__":
    METATLAS_DB_CONFIG = '/home/bkrull/.fireworks/metatlas.ini'
    CSV_FILE = './chebi_and_metacyc_molecules.csv'
    PROJECT_HOME = 'scr/'

    metatlas_lpad = create_launchpad(METATLAS_DB_CONFIG)
    metatlas_lpad.reset('2018-03-13')
    molecules = make_df_with_smiles_only_from_csv(CSV_FILE, reset=True)

    for index, row in molecules.iterrows():
        smiles = row['original_smiles']
        formula = row['formula']

        if type(smiles) is float:
            continue

        uff = OBUFFOptimize(smiles_string=smiles)
        pm3 = OrcaOptimize()
        fw = Firework([uff, pm3], name=formula)

        metatlas_lpad.add_wf(fw)
