from __future__ import annotations

from drafter import *
from dataclasses import dataclass
import inspect


#CHEMICAL REACTION NAVIGATION

def chemical_to_int(chemical : str) -> list[list[int]]:
    ret = []
    i = 0
    while i < len(chemical):
        coeff = 1
        atom = chemical[i]
        if(i < len(chemical) - 1 and chemical[i + 1].isdigit()):
            coeff = chemical[i+1]
        if(atom.isalpha() and atom.upper() == atom):
            if(i < len(chemical) - 1 and chemical[i + 1].isalpha() and chemical[i+1].lower() == chemical[i+1]):
                atom = chemical[i:i+2]
                if(i < len(chemical) - 2):
                    coeff = chemical[i+2]
                else:
                    coeff = 1
        else:
            i += 1
            continue

        number = pt[atom]
        ret.append([number, int(coeff)])
        i += 1
    return ret

def int_to_chemical(chem : list[list[int]]) -> str:
    ret = ""
    for i in chem:
        ret += ptr[i[0]]
        if(i[1] != 1):
            ret += str(i[1])
    return ret

def int_to_reaction(reac : list[list[list[int, list[int]]]]) -> str:
    ret = ""
    for i in reac[0]:
        if(i[0] != 1):
            ret += str(i[0]) + "⋅"
        ret += int_to_chemical(i[1:])
        ret += " + "
    ret = ret[:-2]
    ret += "-> "
    for i in reac[1]:
        if(i[0] != 1):
            ret += str(i[0]) + "⋅"
        ret += int_to_chemical(i[1:])
        ret += " + "
    return ret[:-3]

def reaction_to_int(reac : str) -> list[list[list[int, list[int]]]]:
    ret = [[], []]
    flip = False
    while(len(reac) > 0):
        try:
            chem = reac[:reac.index(" ")]
            reac = reac[reac.index(" ") + 1:]
        except:
            chem = reac
            reac = ""
        if("->" in chem):
            flip = True
            continue
        try:
            coeff = chem[:chem.index("⋅")]
            chem = chem[chem.index("⋅") + 1:]
            int_chem = [int(coeff)] + chemical_to_int(chem)
        except:
            int_chem = [1] + chemical_to_int(chem)

        if not flip and int_chem != [1]:
            ret[0].append(int_chem)
        if flip and int_chem != [1]:
            ret[1].append(int_chem)
    return ret

def equation_balanced(reaction) -> bool:
    if(type(reaction) == type("a")):
        reaction = reaction_to_int(reaction)
    reac_mult = 1
    for i in reaction[0]:
        mult = i[0]
        for j in i[1:]:
            reac_mult *= primes[j[0]] ** (j[1] * mult)

    prod_mult = 1
    for i in reaction[1]:
        mult = i[0]
        for j in i[1:]:
            prod_mult *= primes[j[0]] ** (j[1] * mult)

    return reac_mult == prod_mult

def reaction_hf(reaction) -> float:
    if(type(reaction) == type("a")):
        reaction = reaction_to_int(reaction)
    hfreac = 0.0
    for i in reaction[0]:
        chemical = i[1:]
        chstr = int_to_chemical(chemical)
        if(chstr == None):
            continue
        hfreac -= hf(chstr)

    for i in reaction[1]:
        chemical = i[1:]
        chstr = int_to_chemical(chemical)
        if(chstr == None):
            continue
        hfreac += hf(chstr)
    
    return hfreac


#MATH




#SITE

@dataclass
class Chemical:
    formula: str
    atoms : list[list[int]]
    atoms : list[list[int]]
    heat_of_formation : float

@dataclass
class State:
    reaction_equation : list[list[list[int, list[int]]]]
    display_equation : list[list[int, Chemical]]
    enthalpy : float
    t_thresh : float
    t_room : float
    reactor_volume : float
    flow_rates : list[Chemical, float]

@route
def index(state : State):
    return Page(state, "hello")






#VERIFICATION AND PRINT-TESTING:
test = False

if(test):
    print(hf("H2O"))

    print(chemical_to_int("CH4"))

    print(int_to_chemical([[6, 2], [1, 4]]))

    print(int_to_reaction(
        [
        [ [1, [6, 1], [1, 4]], [2, [8, 2]]],
        [ [1, [6, 1], [8, 2]], [2, [1, 2], [8, 1]]]
        ]
        ))

    print(reaction_to_int("CH4 + 2⋅O2 -> CO2 + 2⋅H2O"))

    print(int_to_reaction(reaction_to_int("CH4 + 2⋅O2 -> CO2 + 2⋅H2O")))

    print(equation_balanced("CH4 + 2⋅O2 -> CO2 + 2⋅H2O"))

    print(equation_balanced(
        [
        [ [1, [6, 1], [1, 4]], [2, [8, 2]]],
        [ [1, [6, 1], [8, 2]], [2, [1, 2], [8, 1]]]
        ]))

    print(reaction_hf("CH4 + 2⋅O2 -> CO2 + 2⋅H2O"))

start_server(State([], [], 0.0, 0.0, 0.0, 0.0, []))




pt = {
    # Period 1
    "H": 1, "He": 2,
    
    # Period 2
    "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10,
    
    # Period 3
    "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15, "S": 16, "Cl": 17, "Ar": 18,
    
    # Period 4
    "K": 19, "Ca": 20, "Sc": 21, "Ti": 22, "V": 23, "Cr": 24, "Mn": 25, "Fe": 26,
    "Co": 27, "Ni": 28, "Cu": 29, "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34,
    "Br": 35, "Kr": 36,
    
    # Period 5
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43, "Ru": 44,
    "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50, "Sb": 51, "Te": 52,
    "I": 53, "Xe": 54,
    
    # Period 6
    "Cs": 55, "Ba": 56, "La": 57, "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62,
    "Eu": 63, "Gd": 64, "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70,
    "Lu": 71, "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85, "Rn": 86,
    
    # Period 7
    "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92, "Np": 93, "Pu": 94,
    "Am": 95, "Cm": 96, "Bk": 97, "Cf": 98, "Es": 99, "Fm": 100, "Md": 101, "No": 102,
    "Lr": 103, "Rf": 104, "Db": 105, "Sg": 106, "Bh": 107, "Hs": 108, "Mt": 109, "Ds": 110,
    "Rg": 111, "Cn": 112, "Nh": 113, "Fl": 114, "Mc": 115, "Lv": 116, "Ts": 117, "Og": 118
}

ptr = [
    None,  # Index 0 placeholder
    
    # Period 1
    "H", "He",
    
    # Period 2
    "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    
    # Period 3
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar",
    
    # Period 4
    "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe",
    "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se",
    "Br", "Kr",
    
    # Period 5
    "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te",
    "I", "Xe",
    
    # Period 6
    "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm",
    "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt",
    "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn",
    
    # Period 7
    "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
    "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
]

primes = [
    2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
    127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
    179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
    233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
    283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
    353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
    419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
    467, 479, 487, 491, 499, 503, 509, 521, 523, 541,
    547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
    607, 613, 617, 619, 631, 641, 643, 647, 653, 659,
    661, 673, 677, 683, 691, 701, 709, 719, 727, 733,
    739, 743, 751, 757, 761, 769, 773, 787, 797, 809,
    811, 821, 823, 827, 829, 839, 853, 857, 859, 863,
    877, 881, 883, 887, 907, 911, 919, 929, 937, 941,
    947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013,
    1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069,
    1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151,
    1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223
]

