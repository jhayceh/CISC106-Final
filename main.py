from drafter import *
from dataclasses import dataclass
import inspect
#CHEMICAL REACTION NAVIGATION


@dataclass
class Chemical:
    formula: str
    atoms : list[list[int]]
    heat_of_formation : float

def chemical_to_int(chemical : str) -> list[list[int]]:
    '''
    Function that takes a string representing a chemical formula and converts
    it into a compound integer representation.

    Params:
        chemical (str): String representing a chemical formula (like CH4)

    Returns:
        compount list of integers representing atoms and quantities in chemical
        (ex CH4 -> [[6, 1], [1, 4]])
    '''
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
    '''
    Function that takes a conpound list of integers representing a chemical formula and converts
    it into a string.

    Params:
        chemical (list[list[int]]): list of ints representing a chemical formula ([[6, 1], [1, 4]])

    Returns:
        string representing chemical
        (ex [[6, 1], [1, 4]] -> CH4)
    '''
    ret = ""
    for i in chem:
        ret += ptr[i[0]]
        if(i[1] != 1):
            ret += str(i[1])
    return ret

def int_to_reaction(reac : list[list[list[int, list[int]]]]) -> str:
    """
    Function that takes a compound list of ints and turns it into a string
    representing a reaction

    Params:
        reac (list[list[list[int, list[int]]]]): reaction as a list

    Returns:
        reaction as a string
        """
    
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
    """
    Function that takes a string and turns it into a compound list of ints
    representing a reaction

    Params:
        reac (str): reaction as a list

    Returns:
        reaction as a list of list of list of int and list of ints representing a reaction
        """
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
    '''
    Takes a reaction and turns it into a list of chemicals because
    I'm required to have a list of dataclass instances.

    Params:
        reac (str): string representing reaction

    Returns:
        List of chemical specices as a list of dataclass instances
    '''
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
    """ 
    Function that checks if a reaction is balanced

    Params:
        reaction (str or list[list[list[int, list[int]]]]): the reaction in
        string or compound list format

    Returns:
        True if the reaction is balanced, else False
    """
    
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





#SITE



@dataclass
class State:
    reaction_equation : str
    involved_species : list[list[Chemical]]
    kin : list[float]
    enthalpy : float
    t_thresh : float
    t_room : float
    reactor_volume : float
    flow_rates : list[Chemical, float]
    flow_q : float = 0.0
    t_max : float = 0.0
    t_min : float = 0.0
    safety_factor : float = 1.0
    rho_cp : float = 4184.0
    ua : float = 0.0
    max_s : float = 0.0

@route
def index(state : State):
    """
    Route creating home page of the site
    """
    return Page(state,
                [
                    Header("Chemical Reactor Enthalpy Calculator"),
                    Button("Start", "builder")
                ])

@route
def builder(state : State):
    """
    Route creating reaction builder page
    """
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
    """
    Route creating route confirmation page
    """
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
    """
    Route creating enthalpy input page
    """
    return Page(state, [
        "Input Enthalpy in kJ per mol reaction: \n\n",
        TextBox("h"),
        "\n",
        Button("Continue", "save_enthalpy")
        ])

@route
def kinetics_data(state : State):
    """
    Route creating kinetics data input page
    """
    rint = reaction_to_int(state.reaction_equation)
    pagec = ["Kinetics Equation Data:\n\nRate = k * ", "\nk:", TextBox("k")]
    for v, i in enumerate(rint[0]):
        crep = (int_to_chemical(i[1:]))
        pagec[0] += "[" + crep + "]^n"
        if(v != len(rint) - 1):
            pagec[0] += " * "
        pagec.append("\nn for " + crep + ":")
        pagec.append(TextBox("che" + str(v)))
    pagec.append("\n\n")
    pagec.append(Button("Go", "kinfin"))                
    return Page(state, pagec)

@route
def kinfin(state : State, k : str, che0 : str = "", che1 : str = "", che2 : str = "", che3 : str = "", che4 : str = "", che5 : str = "", che6 : str = "", che7 : str = ""):
    """
    Route handling kinetics data
    """
    coeff = []
    unsc = [k, che0, che1, che2, che3, che4, che5, che6, che7]
    for i in unsc:
        if(i != ""):
            coeff.append(i)
    fcoeff = []
    for i in coeff:
        fcoeff.append(float(i))
    state.kin = fcoeff
    return reactor_info(state)





#REACTOR MODEL
def reactant_data(state):
    '''
    Function that associates kinetic constants with their chemicals
    '''
    rint = reaction_to_int(state.reaction_equation)
    data = []
    for v, chem in enumerate(rint[0]):
        order = state.kin[v + 1] if v + 1 < len(state.kin) else 1.0
        data.append((int_to_chemical(chem[1:]), chem[0], order))
    return data

def effective_limit(state) -> float:
    '''
    Function that determines which temperature limit to use based on
    reaction enthalpy direction and FoS
    '''
    sf = max(state.safety_factor, 1.0)
    if state.enthalpy < 0:
        return state.t_room + (state.t_max - state.t_room) / sf
    return state.t_room - (state.t_room - state.t_min) / sf

def temp_exceeds(state, T : float) -> bool:
    '''
    Function that checks if the temperature already exceeds the
    threshhold value
    '''
    lim = effective_limit(state)
    if state.enthalpy < 0:
        return T > lim
    return T < lim

def simulate(state, s : float, steps : int = 500, check : bool = True):
    """
    Function that actually handles the simulation of the reactor environment
    """
    species = reactant_data(state)
    k = state.kin[0]
    V = state.reactor_volume
    q = state.flow_q
    rcp = state.rho_cp
    ua = state.ua
    dH = state.enthalpy * 1000.0
    T_amb = state.t_room
    t_end = 10.0 * V / q
    dt = t_end / steps
    cin = [sp[1] * s / q for sp in species]
    n_sp = len(species)

    def deriv(y):
        '''
        subfunction handling the math with the kinetics equation
        '''
        C = y[1:]
        r = k
        for j in range(n_sp):
            c = max(C[j], 0.0)
            n = species[j][2]
            if n < 0 and c <= 0:
                r = 0.0
            elif n != 0:
                r *= c ** n
        dC = []
        for j in range(n_sp):
            dC.append((cin[j] - C[j]) * q / V - species[j][1] * r)
        dT = (q * rcp * (T_amb - y[0]) + (-dH) * r * V - ua * (y[0] - T_amb)) / (rcp * V)
        return [dT] + dC

    y = [T_amb] + [0.0] * n_sp
    times = [0.0]
    temps = [T_amb]
    hit = temp_exceeds(state, T_amb)
    if check and hit:
        return times, temps, True
    size = len(y)
    for step in range(steps):
        k1 = deriv(y)
        k2 = deriv([y[i] + 0.5 * dt * k1[i] for i in range(size)])
        k3 = deriv([y[i] + 0.5 * dt * k2[i] for i in range(size)])
        k4 = deriv([y[i] + dt * k3[i] for i in range(size)])
        y = [y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(size)]
        times.append((step + 1) * dt)
        temps.append(y[0])
        if temp_exceeds(state, y[0]):
            hit = True
            if check:
                return times, temps, True
    return times, temps, hit

def find_max_s(state) -> float:
    """
    Function that finds the largest s (mol/s) that keeps T inside the limit.
    Returns -1.0 if no limit is ever reached.
    """
    if state.enthalpy == 0:
        return -1.0
    if simulate(state, 0.0)[2]:
        return 0.0                        
    lo = 0.0
    hi = 1e-9
    found = False
    for _ in range(60):
        if simulate(state, hi)[2]:
            found = True
            break
        lo = hi
        hi *= 4.0
    if not found:
        return -1.0
    for _ in range(30):                  
        mid = (lo + hi) / 2.0
        if simulate(state, mid)[2]:
            hi = mid
        else:
            lo = mid
    return lo

def text_graph(state, times, temps) -> str:
    """
    Function that creates a text graph of the function
    """
    lim = effective_limit(state)
    lo = min(min(temps), lim, state.t_room)
    hi = max(max(temps), lim, state.t_room)
    span = hi - lo
    if span == 0:
        span = 1.0
    out = "Temperature vs time (bars scaled from " + str(round(lo, 1)) + " to " + str(round(hi, 1)) + " C)\n"
    stride = max(1, len(times) // 25)
    for i in range(0, len(times), stride):
        bar = "#" * int(40 * (temps[i] - lo) / span)
        out += "t=" + str(round(times[i], 1)) + " s  T=" + str(round(temps[i], 2)) + "  " + bar + "\n"
    out += "Effective limit: " + str(round(lim, 2)) + " C"
    return out


#PAGES: REACTOR INFO, RESULTS, GRAPH

@route
def save_enthalpy(state : State, h : str = ""):
    """
    Route called from the continue button on the enthalpy page that saves enthalpy data
    """
    try:
        state.enthalpy = float(h)
    except ValueError:
        return Page(state, ["Enthalpy must be a number (kJ/mol of reaction as written, negative = exothermic).\n",
                            Button("Back", "ent")])
    return kinetics_data(state)

@route
def reactor_info(state : State):
    """
    Route that constructs the page to input reactor and reaction info
    """
    
    return Page(state, [
        "Reactor Information\n\n",
        "Reactor volume (L): ", TextBox("volume"),
        "\nTotal volumetric flow through reactor (L/s): ", TextBox("flow"),
        "\nRoom / feed temperature (C): ", TextBox("troom"),
        "\nMaximum allowed temperature (C): ", TextBox("tmax"),
        "\nMinimum allowed temperature (C): ", TextBox("tmin"),
        "\nFactor of safety (1 or more): ", TextBox("sf"),
        "\nVolumetric heat capacity of contents, rho*Cp (J/L/K) (blank = 4184, water): ", TextBox("rhocp"),
        "\nHeat loss coefficient UA to surroundings (W/K) (blank = 0, adiabatic): ", TextBox("ua"),
        "\n\n", Button("Go to Results", "reactor_save")
    ])

@route
def reactor_save(state : State, volume : str = "", flow : str = "", troom : str = "", tmax : str = "",
                 tmin : str = "", sf : str = "", rhocp : str = "", ua : str = ""):
    """
    Route function that saves reactor data inputted in reactor_info
    """

    try:
        v = float(volume)
        q = float(flow)
        t_room = float(troom)
        t_max = float(tmax)
        t_min = float(tmin)
        safety = float(sf)
        rcp = float(rhocp) if rhocp.strip() != "" else 4184.0
        ua_val = float(ua) if ua.strip() != "" else 0.0
    except ValueError:
        return Page(state, ["Every field must be a number (only the last two may be left blank).\n",
                            Button("Back", "reactor_info")])
    if v <= 0 or q <= 0 or rcp <= 0 or ua_val < 0 or safety < 1:
        return Page(state, ["Volume, flow and rho*Cp must be positive, UA cannot be negative, and the factor of safety must be at least 1.\n",
                            Button("Back", "reactor_info")])
    if t_max <= t_min:
        return Page(state, ["Maximum temperature must be greater than minimum temperature.\n",
                            Button("Back", "reactor_info")])
    state.reactor_volume = v
    state.flow_q = q
    state.t_room = t_room
    state.t_max = t_max
    state.t_min = t_min
    state.safety_factor = safety
    state.rho_cp = rcp
    state.ua = ua_val
    state.t_thresh = t_max if state.enthalpy < 0 else t_min
    return results(state)

@route
def results(state : State):
    """
    Route function that constructs the results page
    """
    if len(state.kin) == 0 or state.reaction_equation == "" or state.reactor_volume <= 0:
        return Page(state, ["Some inputs are missing. Please start over.\n", Button("Restart", "restart")])

    if state.enthalpy == 0:
        return Page(state, ["The reaction has zero enthalpy, so it never heats or cools the reactor. There is no temperature limit on flow rates.\n",
                            Button("Restart", "restart")])

    kind = "exothermic (limit is the MAXIMUM temperature)" if state.enthalpy < 0 else "endothermic (limit is the MINIMUM temperature)"
    state.max_s = find_max_s(state)
    lim = effective_limit(state)

    content = ["Results\n\n",
               "Reaction: " + state.reaction_equation + "\n",
               "Type: " + kind + "\n",
               "Threshold from reactor info: " + str(state.t_thresh) + " C\n",
               "Effective limit with factor of safety " + str(state.safety_factor) + ": " + str(round(lim, 2)) + " C\n\n"]

    if state.max_s < 0:
        content.append("No flow rate in the tested range reaches the limit, so this reactor is not temperature-limited.\n\n")
    elif state.max_s == 0:
        content.append("The limit is already violated with no reaction (room temperature is outside the allowed range after the factor of safety). Flow rates of 0 are the only safe choice.\n\n")
    else:
        content.append("Maximum safe feed rates:\n")
        for formula, nu, order in reactant_data(state):
            content.append(formula + ": " + str(nu * state.max_s) + " mol/s\n")
        times, temps, hit = simulate(state, state.max_s, 500, False)
        if state.enthalpy < 0:
            content.append("\nPeak reactor temperature at these rates: " + str(round(max(temps), 2)) + " C\n\n")
        else:
            content.append("\nLowest reactor temperature at these rates: " + str(round(min(temps), 2)) + " C\n\n")
    content.append(Button("View Graph", "graph"))
    content.append(Button("Restart", "restart"))
    return Page(state, content)

@route
def graph(state : State):
    """
    Route function that constructs the graph page
    """
    
    if len(state.kin) == 0 or state.reactor_volume <= 0 or state.enthalpy == 0:
        return Page(state, ["Nothing to graph yet.\n", Button("Restart", "restart")])
    s = max(state.max_s, 0.0)
    times, temps, hit = simulate(state, s, 500, False)
    lim = effective_limit(state)
    buttons = ["\n", Button("Back to Results", "results"), Button("Restart", "restart")]
    try:
        import matplotlib.pyplot as plt
        plt.figure()
        plt.plot(times, temps, label="Reactor temperature")
        plt.axhline(lim, color="red", linestyle="--", label="Effective limit")
        plt.xlabel("Time (s)")
        plt.ylabel("Temperature (C)")
        plt.title("Reactor temperature at maximum safe flow")
        plt.legend()
        return Page(state, ["Graph of reactor temperature over time (flow rates from Results).\n", MatPlotLibPlot()] + buttons)
    except Exception:
        return Page(state, ["Graph of reactor temperature over time (flow rates from Results).\n",
                            text_graph(state, times, temps)] + buttons)

@route
def restart(state : State):
    """
    Function that returns state to what it was at the start of the program and
    reruns the program
    """
    state.reaction_equation = ""
    state.involved_species = []
    state.kin = []
    state.enthalpy = 0.0
    state.t_thresh = 0.0
    state.t_room = 0.0
    state.reactor_volume = 0.0
    state.flow_rates = []
    state.flow_q = 0.0
    state.t_max = 0.0
    state.t_min = 0.0
    state.safety_factor = 1.0
    state.rho_cp = 4184.0
    state.ua = 0.0
    state.max_s = 0.0
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


#PERIODIC TABLE AND PRIMES

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



start_server(State("", [], [], 0.0, 0.0, 0.0, 0.0, []))
