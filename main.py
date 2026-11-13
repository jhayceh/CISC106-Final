from __future__ import annotations

from drafter import *
from dataclasses import dataclass
import inspect
from parse import pt, ptr, primes


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

def reaction_to_chemicals(reac : str) -> list[list[int, Chemical]]:
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

        add_chem = Chemical(chem, int_chem, 0.0)

        if not flip and int_chem != [1]:
            ret[0].append(add_chem)
        if flip and int_chem != [1]:
            ret[1].append(add_chem)
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


#MATH




#SITE

@dataclass
class Chemical:
    formula: str
    atoms : list[list[int]]
    heat_of_formation : float

@dataclass
class State:
    reaction_equation : str
    involved_species : list[list[Chemical]]
    enthalpy : float
    t_thresh : float
    t_room : float
    reactor_volume : float
    flow_rates : list[Chemical, float]

@route
def index(state : State):
    return Page(state,
                [
                    Header("Chemical Reactor Enthalpy Calculator"),
                    Button("Start", "builder")
                ])

@route
def builder(state : State):
    return Page(state,
                [
                    "Input your reaction:\n\n",
                    "Reactants: ", TextBox("left"),
                    "\nProducts: ", TextBox("right"),
                    "\n", Button("Confirm", "creac")
                    ]
                )

@route
def creac(state : State, left, right):
    reaction = left + " -> " + right
    i = 0;
    while i < len(reaction):
        if(reaction[i].isdigit() and reaction[i+1].isalpha() and (i == 0 or reaction[i-1].isdigit() or reaction[i-1] == " ")):
            reaction = reaction[:i+1] + "⋅" + reaction[i+1:]
            i += 1
        i += 1
    rint = reaction_to_int(reaction)
    state.reaction_equation = reaction
    return Page(state,
                [
                    "Please confirm this is your recation:\n",
                    reaction,
                    "\n",
                    Button("Back", "builder"), Button("Continue", "ent")
                    ]
                )

@route
def ent(state : State):
    return Page(state, [
        "Input Enthalpy: \n\n",
        TextBox("h"),
        "\n",
        Button("Continue", "data_in")
        ])

@route
def data_in(state : State):
    return index(state)

    




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

start_server(State("", [], 0.0, 0.0, 0.0, 0.0, []))
