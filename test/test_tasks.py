#! /usr/bin/env python
import pytest
from metatlas import create_pybel_molecule, OBUFFOptimize


@pytest.mark.parametrize('mol_str, strtype, formula', [
    ('O[C@H]1C[C@@]23C[C@]2(C1)OC3', 'smi', 'C7H10O2'),
    ('/home/bkrull/Documents/data/qm9/dsgdb9nsd_090177.xyz', '', 'C7H10O2'),
    ('/home/bkrull/Documents/data/qm9/dsgdb9nsd_090177.xyz', 'xyz', 'C7H10O2')
])
def test_make_pybel_molecule(mol_str, strtype, formula):
    if strtype == '':
        mol = create_pybel_molecule(mol_str)
    else:
        mol = create_pybel_molecule(mol_str, strtype)

    if mol.formula == formula:
        assert True
    else:
        assert False, "Test formula: {}; expected {}".format(mol.formula,
                                                             formula)


specs = [({
    '_fw_name': 'OBUFFOptimize', 'strtype': 'xyz', '_expected': 'C7H10O2',
    'smiles_string': '/home/bkrull/Documents/data/qm9/dsgdb9nsd_090177.xyz',
}),({
    '_fw_name': 'OBUFFOptimize', 'strtype': 'smi', '_expected': 'C7H10O2',
    'smiles_string': 'O[C@H]1C[C@@]23C[C@]2(C1)OC3',
}),({
    '_fw_name': 'OBUFFOptimize', 'strtype': '', '_expected': 'C7H10O2',
    'smiles_string': '/home/bkrull/Documents/data/qm9/dsgdb9nsd_090177.xyz',
})]
@pytest.mark.parametrize('spec', specs)
def test_obuff(spec):

    fwaction = OBUFFOptimize(spec)
    result = fwaction.run_task(spec).to_dict()

    assert result['mod_spec'][0]['_push']['formula'] == \
        spec['_expected'], 'Result value: {}; expected: {}'
