import sys
from generation_steps import *
sys.path.append("..")


def write_out(tmp, name):
    """
    write smiles in set to file
    :param tmp:
    :param name:
    :return:
    """
    with open(name, 'w') as out:
        for smi in tmp:
            out.write(smi + "\n")


def func_dict(step_name):
    my_dict = {
        'step_1': add_ring,
        'step_2': sidechain_walk,
        'step_3': add_aromatic_ring,
        'step_4': bond_walk,
        'step_5': hetero_walk,
        'step_6': others_walk,
        'step_7': bond_walk,
        'step_8': sidechain_walk
    }

    return my_dict[step_name]


def process(num):
    """
    generate all possible smiles with fixed heavy atom number
    :param num:
    :return:
    """
    target = 'C' * num
    out = {target}
    for i in range(1, 8):
        func = func_dict('step_' + str(i))
        out = func(out)
    filename = str(num) + "_out.txt"
    write_out(out, filename)


def main():
    for i in range(3, 13):
        process(i)
