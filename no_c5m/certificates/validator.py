import pickle
from math import comb, factorial
from fractions import Fraction
from itertools import permutations, combinations
import argparse


class flag:
    def __init__(self, size, edges=[], ftype=[], parts=[], col2=[], col3=[]):
        self.size = size
        self.edges = edges
        self.ftype = ftype
        self.parts = parts
        self.col2 = set(col2)
        self.col3 = set(col3)

    def subflag(self, vertices):
        new_edges = [edge for edge in self.edges
                     if all(v in vertices for v in edge)]
        new_parts = [part for part in self.parts
                     if all(v in vertices for v in part)]
        new_col2 = self.col2 & set(vertices)
        new_col3 = self.col3 & set(vertices)
        relabel = {vertices[i]: i for i in range(len(vertices))}
        new_edges_1 = [[relabel[v] for v in edge] for edge in new_edges]
        new_ftype_1 = [relabel[v] for v in self.ftype]
        new_parts_1 = [[relabel[v] for v in part] for part in new_parts]
        new_col2_1 = [relabel[v] for v in new_col2]
        new_col3_1 = [relabel[v] for v in new_col3]
        return flag(
            size=len(vertices),
            edges=new_edges_1,
            ftype=new_ftype_1,
            parts=new_parts_1,
            col2=new_col2_1,
            col3=new_col3_1
        )


# ---------------------------- Auxiliary functions ----------------------------

def is_isomorphism(bijection, flag1, flag2):
    """Checks if a vertex bijection between two flags is an isomorphism.

    """
    # Check type mapping
    if flag2.ftype != [bijection[v] for v in flag1.ftype]:
        return False
    # Check edge mapping
    mapped_edges = set()
    for edge in flag1.edges:
        mapped_edge = frozenset(bijection[v] for v in edge)
        mapped_edges.add(mapped_edge)
    if mapped_edges != {frozenset(edge) for edge in flag2.edges}:
        return False
    # Check part mapping
    mapped_parts = set()
    for part in flag1.parts:
        mapped_part = frozenset(bijection[v] for v in part)
        mapped_parts.add(mapped_part)
    if mapped_parts != {frozenset(part) for part in flag2.parts}:
        return False
    # Check colour mapping
    if flag2.col2 != {bijection[v] for v in flag1.col2}:
        return False
    if flag2.col3 != {bijection[v] for v in flag1.col3}:
        return False
    return True


def are_isomorphic(flag1, flag2):
    """Checks if two flags are isomorphic.

    """
    # Check number of edges and parts
    if len(flag1.edges) != len(flag2.edges):
        return False
    if len(flag1.parts) != len(flag2.parts):
        return False
    # Check colour classes
    if len(flag1.col2) != len(flag2.col2):
        return False
    if len(flag1.col3) != len(flag2.col3):
        return False
    V = range(flag1.size)
    # Check degree sequence
    if {len([e for e in flag1.edges if v in e]) for v in V} != {len([e for e in flag2.edges if v in e]) for v in V}:
        return False
    # Check every permutation
    for permutation in permutations(V):
        bijection = {v: u for v, u in zip(V, permutation)}
        if is_isomorphism(bijection, flag1, flag2):
            return True
    return False


def change_base(quantum, base_flags):
    """Given a qunatum flag and some base flags of the same size, expresses the
    quantum flag as a linear combination of the base flags. Returns the coefficients
    of this vector.

    """
    coefficients = [Fraction(0) for _ in range(len(base_flags))]
    for coef, small in quantum:
        for idx, big in enumerate(base_flags):
            # Obtain the density of 'small' in 'big'
            copies = 0
            # Check every small-size subset of big for isomorphism
            for U in combinations(range(big.size), small.size):
                if are_isomorphic(big.subflag(U), small):
                    copies += 1
            density = Fraction(copies, comb(big.size, small.size))
            coefficients[idx] += coef * density
    return coefficients


def get_densities(quantum_flag, flags, G):
    """Given a typed quantum flag, a list of typed flags, and a host 0-flag,
    calculates the list of densities of the type-average of the product of the
    quantum flag with each of the flags in the host flag.

    """
    n = G.size
    n_flags = len(flags)
    n_qflag = len(quantum_flag)
    flags_size = flags[0].size
    quantum_size = quantum_flag[0][1].size
    type_size = len(flags[0].ftype)
    V = range(n)
    M = [[0 for _ in range(n_flags)] for _ in range(n_qflag)]
    for ftype in permutations(V, type_size):
        new_flag = flag(G.size, G.edges, list(ftype), G.parts, G.col2, G.col3)
        unlabeled_vertices = [v for v in V if v not in ftype]
        for V1 in combinations(unlabeled_vertices, quantum_size - type_size):
            flag1 = new_flag.subflag(ftype + V1)
            for i, flag_a in enumerate(quantum_flag):
                if are_isomorphic(flag1, flag_a[1]):
                    remaining_vertices = [v for v in unlabeled_vertices if v not in V1]
                    for V2 in combinations(remaining_vertices, flags_size - type_size):
                        flag2 = new_flag.subflag(ftype + V2)
                        for j, flag_b in enumerate(flags):
                            if are_isomorphic(flag2, flag_b):
                                M[i][j] += 1
                                break
                    break
    denominator = comb(n, type_size) * factorial(type_size)
    denominator *= comb(n - type_size, quantum_size - type_size)
    denominator *= comb(n - quantum_size, flags_size - type_size)
    densities = [0 for _ in range(n_flags)]
    for i in range(n_flags):
        numerator = sum([quantum_flag[j][0] * M[j][i] for j in range(n_qflag)])
        densities[i] = Fraction(numerator, denominator)
    return densities


def average(typed_flag):
    """Given a typed flag, calculates the probability that the subjacent 0-flag,
    together with a uniformly at random type is isomorphic to the flag.

    """
    n, edges, parts, col2, col3 = typed_flag.size, typed_flag.edges, typed_flag.parts, typed_flag.col2, typed_flag.col3
    type_size = len(typed_flag.ftype)
    numerator = 0
    for ftype in permutations(range(n), type_size):
        if are_isomorphic(typed_flag, flag(n, edges, list(ftype), parts, col2, col3)):
            numerator += 1
    return Fraction(numerator, comb(n, type_size) * factorial(type_size))


def recover_from_certificate(file_name):
    """Given a certificate, recovers the stored objects.

    """
    with open(file_name, 'rb') as handle:
        certificate = pickle.load(handle)
        base_flags_0 = certificate['base flags']
        target = certificate['target']
        positives = certificate['positives']
        factor_flags_0 = certificate['factor flags']

    # Get list of base flags and typed flags in format
    base_flags = [flag(size=f[0],
                       edges=[list(e) for e in f[2][0][1]], ftype=[],
                       parts=[list(x) for x in (f[2][1][1]
                                                if len(f[2]) == 2 else [])],
                       col2=[x[0] for x in (f[2][2][1] if len(f[2]) == 4 else [])],
                       col3=[x[0] for x in (f[2][3][1] if len(f[2]) == 4 else [])]
                       )
                  for f in base_flags_0]

    factor_flags = []
    for a in factor_flags_0:
        b = []
        for c in a:
            b += [flag(c[0],
                       [list(e) for e in c[2][0][1]],
                       list(c[1]),
                       [list(x) for x in (c[2][1][1]
                                                   if len(c[2]) == 2 else [])],
                       [x[0] for x in (c[2][2][1] if len(c[2]) == 4 else [])],
                       [x[0] for x in (c[2][3][1] if len(c[2]) == 4 else [])]
                       )
                  ]
        factor_flags += [b]
    return base_flags, target, factor_flags, positives


# ---------------------------- Objectives and assumptions ----------------------------

def prop_3_1():
    """Returns the objective function as a quantum flag and the positivity
    assumptions as a list of quantum flags of Proposition 3.1.

    """
    objective = [(1, flag(3, [[0, 1, 2]]))]
    positivity_assumptions = []
    return objective, positivity_assumptions


def prop_3_2(base_flags):
    """Returns the objective function as a quantum flag and the positivity
    assumptions as a list of quantum flags of Proposition 3.2.

    """
    # Define objective function
    print("Defining objective function...")

    # Get list of all 6-vertex flags from the 7-vertex base flags
    six_flags = []
    for f in base_flags:
        for v in range(7):
            subflag = f.subflag([x for x in range(7) if x != v])
            if all(not are_isomorphic(r, subflag) for r in six_flags):
                six_flags += [subflag]

    # Define the typed quantum flag
    base = [[0, 1, 2], [3, 4, 5]]
    f222 = [
        flag(6, edges=base+[[0, 1, 3], [0, 2, 5], [1, 2, 4]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 1, 3], [0, 3, 5], [0, 2, 5], [1, 2, 4]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 1, 3], [0, 3, 5], [0, 2, 5], [1, 2, 4]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 1, 3], [0, 3, 5], [0, 2, 5], [1, 2, 4]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 1, 3], [0, 3, 5], [0, 2, 5], [1, 2, 4], [2, 4, 5]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 2, 3], [0, 3, 5], [0, 1, 5], [1, 2, 4], [1, 4, 5]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 4, 5], [2, 3, 4], [0, 2, 4], [0, 1, 5], [1, 2, 3]], ftype=[0, 1, 2]),
        flag(6, edges=base+[[0, 2, 4], [0, 1, 5], [0, 4, 5], [1, 2, 3], [2, 3, 4], [1, 3, 5]], ftype=[0, 1, 2])
    ]

    # Average the terms to produce a 0-quantum flag using the 6-vertex flags
    objective = []
    for i, G in enumerate(six_flags):
        coefficient = sum([average(term) for term in f222
                           if are_isomorphic(G, flag(term.size, term.edges))])
        if coefficient > 0:
            objective += [(coefficient, G)]

    # Lower bound on average degree assumption (edge - beta_3_2 >= 0)
    beta_3_2 = Fraction(1, 4) - Fraction(1, 100000)
    avg_degree_minus_beta_3_2 = [
        (-beta_3_2, flag(3)),
        (1-beta_3_2, flag(3, [[0, 1, 2]]))
    ]
    positivity_assumptions = [avg_degree_minus_beta_3_2]
    return objective, positivity_assumptions


def prop_3_3():
    """Returns the objective function as a quantum flag and the positivity
    assumptions as a list of quantum flags of Proposition 3.3.

    """
    objective = [
        (1, flag(3, edges=[[0, 1, 2]], parts=[[0, 1], [1, 2]])),
        (-Fraction(99, 100), flag(3, parts=[[0, 1], [0, 2], [1, 2]]))
    ]

    # Define assumptions

    # Local optimality
    Cp_minus_Bp_2 = [
        (1, flag(3, edges=[[0, 1, 2]], ftype=[0], parts=[[0, 1], [0, 2], [1, 2]])),
        (-Fraction(1, 2), flag(3, edges=[[0, 1, 2]], ftype=[0], parts=[[0, 1], [1, 2]]))
    ]

    # Edge density is larger than beta_3_3
    beta_3_3 = Fraction(19, 100)
    C_minus_beta_3_3 = [
        (-beta_3_3, flag(3)),
        (-beta_3_3, flag(3, parts=[[0, 1], [0, 2]])),
        (-beta_3_3, flag(3, parts=[[0, 1], [0, 2], [1, 2]])),
        (-beta_3_3, flag(3, edges=[[0, 1, 2]])),
        (-beta_3_3, flag(3, edges=[[0, 1, 2]], parts=[[0, 1], [0, 2]])),
        (1-beta_3_3, flag(3, edges=[[0, 1, 2]], parts=[[0, 1], [0, 2], [1, 2]]))
    ]

    positivity_assumptions = [Cp_minus_Bp_2, C_minus_beta_3_3]
    return objective, positivity_assumptions


def prop_3_4():
    """Returns the objective function as a quantum flag and the positivity
    assumptions as a list of quantum flags of Proposition 3.4.

    """
    objective = [
        (1, flag(2, edges=[[0, 1]], col2=[1], col3=[])),
        (1, flag(2, edges=[[0, 1]], col2=[], col3=[1])),
        (1, flag(2, edges=[[0, 1]], col2=[0, 1], col3=[])),
        (1, flag(2, edges=[[0, 1]], col2=[], col3=[0, 1])),
        (-Fraction(9, 10), flag(2, col2=[0], col3=[1]))
    ]

    edge_12_minus_01 = [
        (1, flag(2, edges=[[0, 1]], col2=[0], col3=[1])),
        (-1, flag(2, edges=[[0, 1]], col2=[1]))
    ]

    edge_12_minus_02 = [
        (1, flag(2, edges=[[0, 1]], col2=[0], col3=[1])),
        (-1, flag(2, edges=[[0, 1]], col3=[1]))
    ]

    edge_01_02_12_minus_1_8 = [
        (-Fraction(1, 8), flag(2, edges=[], col2=[], col3=[0, 1])),
        (-Fraction(1, 8), flag(2, edges=[], col2=[0], col3=[1])),
        (-Fraction(1, 8), flag(2, edges=[], col2=[0, 1], col3=[])),
        (-Fraction(1, 8), flag(2, edges=[], col2=[], col3=[1])),
        (-Fraction(1, 8), flag(2, edges=[], col2=[1], col3=[])),
        (-Fraction(1, 8), flag(2, edges=[], col2=[], col3=[])),
        (-Fraction(1, 8), flag(2, edges=[[0, 1]], col2=[], col3=[0, 1])),
        (Fraction(7, 8), flag(2, edges=[[0, 1]], col2=[0], col3=[1])),
        (-Fraction(1, 8), flag(2, edges=[[0, 1]], col2=[0, 1], col3=[])),
        (Fraction(7, 8), flag(2, edges=[[0, 1]], col2=[], col3=[1])),
        (Fraction(7, 8), flag(2, edges=[[0, 1]], col2=[1], col3=[])),
        (-Fraction(1, 8), flag(2, edges=[[0, 1]], col2=[], col3=[])),
    ]

    point_0_minus_1_4 = [
        (-Fraction(1, 4), flag(1, col3=[0])),
        (-Fraction(1, 4), flag(1, col2=[0])),
        (Fraction(3, 4), flag(1))
    ]

    point_1_minus_1_4 = [
        (-Fraction(1, 4), flag(1, col3=[0])),
        (Fraction(3, 4), flag(1, col2=[0])),
        (-Fraction(1, 4), flag(1))
    ]

    point_2_minus_1_4 = [
        (Fraction(3, 4), flag(1, col3=[0])),
        (-Fraction(1, 4), flag(1, col2=[0])),
        (-Fraction(1, 4), flag(1))
    ]

    positivity_assumptions = [
        edge_12_minus_01,
        edge_12_minus_02,
        edge_01_02_12_minus_1_8,
        point_0_minus_1_4,
        point_1_minus_1_4,
        point_2_minus_1_4
    ]
    return objective, positivity_assumptions


# ---------------------------- Main ----------------------------

# Read certificate
# =================================================
parser = argparse.ArgumentParser(
    description="Validates that the objectives and assumptions of a given certificate correspond to those in the corresponding proposition.",
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="""
examples:
  $ python validator.py proposition_3_1.pickle
""")
parser.add_argument("certificate", type=str, help='path of the pickle file',
                    choices=['proposition_3_1.pickle', 'proposition_3_2.pickle',
                             'proposition_3_3.pickle', 'proposition_3_4.pickle']
                    )

certificate_path = parser.parse_args().certificate

(base_flags,
 target,
 factor_flags,
 positives) = recover_from_certificate(certificate_path)

# Obtain objective and assumptions
# =================================================
if certificate_path == 'proposition_3_1.pickle':
    objective, positivity_assumptions = prop_3_1()
elif certificate_path == 'proposition_3_2.pickle':
    objective, positivity_assumptions = prop_3_2(base_flags)
elif certificate_path == 'proposition_3_3.pickle':
    objective, positivity_assumptions = prop_3_3()
elif certificate_path == 'proposition_3_4.pickle':
    objective, positivity_assumptions = prop_3_4()


# Validate objective
# =================================================
if isinstance(objective[0], tuple):
    objective = change_base(objective, base_flags)
if target == objective:
    print("Objective function matches.\n")
else:
    print("ERROR: OBJECTIVE FUNCTION DOES NOT MATCH!")
    exit(1)


# Validate assumptions
# =================================================
for idx, assumption in enumerate(positivity_assumptions):
    densities = [get_densities(assumption, factor_flags[idx], G)
                 for G in base_flags]
    # print(len(densities), len(densities[0]))
    transposed = list(map(list, zip(*densities)))
    n_prev_blocks = sum([len(factor_flags[idx]) for idx in range(idx)])
    start, end = n_prev_blocks, n_prev_blocks + len(factor_flags[idx])
    if transposed == positives[start:end]:
        print("Positivity assumption %d matches." % (idx + 1))
    else:
        print("ERROR: POSITIVITY ASSUMPTION %d DOES NOT MATCH!" % (idx + 1))
        exit(1)

print("\nCOMPLETED: THE OBJECTIVES AND ASSUMPTIONS CORRESPOND TO THOSE IN THE PROPOSITION!")
