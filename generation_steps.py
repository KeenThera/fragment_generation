from rdkit import Chem
from copy import deepcopy
from rdkit.Chem.rdChemReactions import *
from rdkit.Chem.EnumerateHeterocycles import EnumerateHeterocycles


def carbon_chain():
    """
    generate sequential carbon string with 1 to 12 'C', such as 'CCCCC'.
    :return:
    """
    out_string = set()
    for leg in range(1, 13):
        tmp = "C" * leg
        out_string.add(tmp)
    return out_string


def gen_index_block(num, block_len):
    """

    :param num:
    :param block_len:
    :return:
    """
    out = []
    for i in range(num):
        if i + block_len <= num:
            out.append((i, i + block_len))
    return out


def get_bridged_atoms(mol):
    """
    collect the bridged atoms in a molecule
    :param mol:
    :return: atom numbers of bridged atoms
    """
    ri = mol.GetRingInfo()
    bond_rings = ri.BondRings()
    bridged_atoms = set()
    for i in range(len(bond_rings)):
        bond_ring_i = set(bond_rings[i])
        for j in range(i):
            bond_ring_j = set(bond_rings[j])
            common_bonds = bond_ring_i.intersection(bond_ring_j)
            if len(common_bonds) > 1:
                atoms = [0 for x in range(len(mol.GetAtoms()))]
                bridged_unit = ()
                for b in common_bonds:
                    atoms[mol.GetBondWithIdx(b).GetBeginAtomIdx()] += 1
                    atoms[mol.GetBondWithIdx(b).GetEndAtomIdx()] += 1
                for idx in range(len(atoms)):
                    if atoms[idx] == 1:
                        bridged_unit += (idx,)
                bridged_atoms.add(bridged_unit)
    return bridged_atoms


def bridged_atoms_filter(mol):
    """
    remove the molecule with aromatic bridged atoms
    :param mol:
    :return:
    """
    bridged_atoms = get_bridged_atoms(mol)
    for i in bridged_atoms:
        for j in i:
            atom = mol.GetAtomWithIdx(j)
            if atom.GetIsAromatic():
                return False
    return True


def substructure_filter(mol):
    """
    remove molecules with unwanted substructures
    :param mol:
    :return:
    """
    # remove ugly ring
    sma_list = ['*1[*R2][*R2]1', '*1[*R2][*R2]*1', '*1[*R2][*R2]*1', '[*x4][r3]']
    for sma in sma_list:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(sma)):
            return False
    # remove molecules with 2 or more [S] atoms
    if len(mol.GetSubstructMatches(Chem.MolFromSmarts('[S]'))) > 1:
        return False

    return True


def max_ring_size_filter(mol):
    """
    remove molecules with max ring size larger than 7
    :param mol:
    :return:
    """
    ring_list = mol.GetRingInfo().AtomRings()
    if len(ring_list) > 0:
        max_ring_size = max([len(ring) for ring in ring_list])
        if max_ring_size > 7:
            return False

    return True


def other_filter(mol):
    """
    remove molecules with heavy atoms larger than 12
    :param mol:
    :return:
    """
    if mol.GetNumAtoms() > 12:
        return False

    return True


def check_mol(smi):
    """
    perform filtering on a given smiles
    :param smi:
    :return:
    """
    mol = Chem.MolFromSmiles(smi)
    if mol is None:
        return False
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False

    if not bridged_atoms_filter(mol):
        return False
    if not substructure_filter(mol):
        return False
    if not other_filter(mol):
        return False
    if not max_ring_size_filter(mol):
        return False

    return True


def ring_generator(smi, depth):
    """
    rings enumeration based on smiles strings
    :param smi:
    :param depth: label numbers (C1CC1), namely ring numbers
    :return:
    """
    start = [(smi, 1)]
    out = set()
    while start:
        cur_smi, cur_dep = start.pop(0)
        if depth is not None and cur_dep > depth:
            return out
        cur_smi = list(cur_smi)
        for block_len in range(2, 2 + len(cur_smi) + 2 * cur_dep):
            for block in gen_index_block(len(cur_smi), block_len):
                new = deepcopy(cur_smi)
                new.insert(block[0], str(cur_dep))
                new.insert(block[1] + 1, str(cur_dep))
                tmp = "".join(new)
                if not check_mol(tmp):
                    continue
                if tmp in out:
                    continue
                start.append((tmp, cur_dep + 1))
                out.add(tmp)


def add_ring(tmp):
    """
    explore new ring based on given smiles set
    :param tmp: smiles set
    :return:
    """
    out = set()
    for smi in tmp:
        if len(smi) < 3:
            continue
        tmp_set = ring_generator(smi, 2)
        if tmp_set is None:
            continue
        out = out.union(tmp_set)
    out = out.union(tmp)
    return out


def my_reaction(rxn_sma, reactants):
    """
    perform reaction on molecules
    :param rxn_sma:
    :param reactants:
    :return:
    """
    rxn = ReactionFromSmarts(rxn_sma)
    ps = rxn.RunReactants(reactants)
    tmp = set([Chem.MolToSmiles(x[0], isomericSmiles=True) for x in ps])
    out = set()
    for i in tmp:
        if Chem.MolFromSmiles(i):
            out.add(i)
    return out


def enumerate_modification(smi, depth, reaction_list):
    mol = Chem.MolFromSmiles(smi)
    input_s = [(mol, 0)]
    seen = set()

    while input_s:
        curmol, dep = input_s.pop(0)
        if depth is not None and dep >= depth:
            return seen
        for rxn_sma in reaction_list:
            rxn = ReactionFromSmarts(rxn_sma)
            for newmol in rxn.RunReactants((curmol,)):
                newmol = newmol[0]
                try:
                    Chem.SanitizeMol(newmol)
                except Exception as e:
                    continue
                newmol_smiles = Chem.MolToSmiles(newmol)
                if newmol_smiles in seen:
                    continue
                input_s.append((newmol, dep + 1))
                seen.add(newmol_smiles)
    return seen


def aromatic_ring_generation_smarts_dictionary():
    """
    smarts dictionary for aromatic ring conversion
    :return:
    """
    runs_dict = {'[C:1]1[C:2][C:3][C:4][C:5][C:6]1>>[c:1]1[c:2][c:3][c:4][c:5][c:6]1': '[*r6]'}
    for i in ['o', 'nH', 's']:
        rxn_sma = "[C:1]1[C:2][C:3][C:4][CH2:5]1>>[c:1]1[c:2][c:3][c:4][{0}:5]1".format(i)
        runs_dict[rxn_sma] = '[*r5]'
        rxn_sma_two = '[C:3]1[C:2]2[C:1]([C:6][C:5][C:4]1)[C:7][C:8][CH2:9]2>>[c:3]1[c:2]2[c:1]' \
                      '([c:6][c:5][c:4]1)[c:7][c:8][{0}:9]2'.format(i)
        runs_dict[rxn_sma_two] = '[CR2r5][CR2r5]'

    runs_dict['[C:1]1[C:2][C:3][C:4][CH1:5]1>>[c:1]1[c:2][c:3][c:4][n:5]1'] = '[*r5]'
    rxn_sma_two_n = '[C:3]1[C:2]2[C:1]([C:6][C:5][C:4]1)[C:7][C:8][CH:9]2>>[c:3]1[c:2]2[c:1]' \
                    '([c:6][c:5][c:4]1)[c:7][c:8][n:9]2'
    runs_dict[rxn_sma_two_n] = '[CR2r5][CR2r5]'
    rxn_sma = '[C:3]1[C:2]2[C:1]([C:6][C:5][C:4]1)[C:7][C:8][C:9][C:10]2>>[c:3]1[c:2]2[c:1]' \
              '([c:6][c:5][c:4]1)[c:7][c:8][c:9][c:10]2'
    runs_dict[rxn_sma] = '[CR2r6][CR2r6]'

    return runs_dict


def aromatic_ring(smi, rxn_dict):
    """
    perform aromatic ring conversion
    :param smi:
    :param rxn_dict:
    :return:
    """
    return enumerate_modification(smi, 2, rxn_dict)


def add_aromatic_ring(tmp):
    """
    explore new aromatic ring on given smiles set
    :param tmp:
    :return:
    """
    out = set()
    my_dict = aromatic_ring_generation_smarts_dictionary()
    for i in tmp:
        new = aromatic_ring(i, my_dict)
        if new is None or len(new) < 1:
            continue
        out = out.union(new)
    out = out.union(tmp)
    return out


def bond_walk(tmp):
    """
    bond transformations on smiles set
    :param tmp:
    :return:
    """
    out = set()
    reaction = ['[A!H0!$(A=*):1][A!H0!$(A=*):2]>>[A:1]=[A:2]', '[A!H0!H1R0:1][A!H0!H1R0:2]>>[A:1]#[A:2]']
    for smi in tmp:
        smi_tmp = enumerate_modification(smi, 1, reaction)
        out = out.union(smi_tmp)

    out = out.union(tmp)
    return out


def hetero_walk(tmp):
    """
    aliphatic hetero atom transformations on given smiles set
    :param tmp:
    :return:
    """
    reaction = ['[C!H0!$(CN):1]>>[N:1]', '[C!H0!H1!$(CO):1]>>[O:1]', '[C!H0!H1!$(CS):1]>>[S:1]']
    out = set()
    for smi in tmp:
        mol = Chem.MolFromSmiles(smi)
        new = enumerate_modification(smi, 3, reaction)
        if new is None:
            print(smi)
            continue
        out = out.union(new)
        for m in EnumerateHeterocycles(mol):
            out.add(Chem.MolToSmiles(m))

    out = out.union(tmp)
    return out


def sidechain_walk(tmp):
    """
    sidechain transformations on given smiles set
    :param tmp:
    :return:
    """
    out = set()
    reaction = ['[A!H0:1][A:2][A:3]>>[A:1]([A:3])[A:2]']
    for smi in tmp:
        new = enumerate_modification(smi, 1, reaction)
        out = out.union(new)

    out = out.union(tmp)
    return out


def others_walk(tmp):
    """
    other transformations on given smiles set
    new valid transformations are welcomed.
    :param tmp: smiles set
    :return:
    """
    out = set()
    reaction = ['[CH3:1]>>[F:1]', '[a:1][CH3,NH2:2]>>[a:1][F:2]']
    for smi in tmp:
        new = enumerate_modification(smi, 1, reaction)
        out = out.union(new)
    out = out.union(tmp)
    return out
