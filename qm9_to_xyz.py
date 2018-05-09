import os
from metatlas import create_launchpad, make_xyz_from_stored_data, make_grad_from_stored_data


if __name__ == '__main__':
    LOCAL_DB_CONFIG = '/home/bkrull/.fireworks/local_db.ini'
    lpad = create_launchpad(LOCAL_DB_CONFIG)

    fract = 1

    for id in lpad.get_fw_ids({'state': 'COMPLETED'}):
        data = lpad.get_fw_dict_by_id(id)
        name = data['name'].split('/')[-1]
        stored_data = data['launches'][0]['action']['stored_data']

        conformers = stored_data['atom_list']
        nconformers = len(conformers)

        try:
            atoms = conformers[0]
        except IndexError:
            continue

        for confid in range(nconformers):
            coords  = stored_data['coords'][confid]
            grads = stored_data['grads'][confid]
            energies = stored_data['energies'][confid]
            with open('/scratch/users/bkrull/qm9struct/{}/{}-{}.xyz'.format(fract, name, confid), 'w') as f:
                for step in range(len(coords)):
                    coord = coords[step]
                    energy = energies[step]
                    xyz = make_xyz_from_stored_data(atoms, coord, energy)
                    f.write(xyz)

            with open('/scratch/users/bkrull/qm9struct/{}/{}-{}.grad'.format(fract, name, confid), 'w') as f:
                for step in range(len(grads)):
                    grad = grads[step]
                    xyz = make_grad_from_stored_data(atoms, grad)
                    f.write(xyz)

        if id % 2500 == 0:
            fract += 1
            os.mkdir('/scratch/users/bkrull/qm9struct/{}'.format(fract))
