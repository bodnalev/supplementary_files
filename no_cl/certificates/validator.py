import pickle
from math import comb, factorial
from fractions import Fraction
from itertools import permutations, combinations
import argparse


class flag:
    def __init__(self, size, edges=[], ftype=[], coloured=[]):
        self.size = size
        self.edges = edges
        self.ftype = ftype
        self.coloured = set(coloured)

    def subflag(self, vertices):
        new_edges = [edge for edge in self.edges
                     if all(v in vertices for v in edge)]
        new_coloured = self.coloured & set(vertices)
        relabel = {vertices[i]: i for i in range(len(vertices))}
        new_edges_1 = [[relabel[v] for v in edge] for edge in new_edges]
        new_ftype_1 = [relabel[v] for v in self.ftype]
        new_coloured_1 = [relabel[v] for v in new_coloured]
        return flag(
            size=len(vertices),
            edges=new_edges_1,
            ftype=new_ftype_1,
            coloured=new_coloured_1
        )


# ---------------------------- Auxiliary functions ----------------------------

def is_isomorphism(bijection, flag1, flag2):
    """Checks if a vertex bijection between two flags is an isomorphism.

    """
    # Check type and colour mapping
    if flag2.ftype != [bijection[v] for v in flag1.ftype]:
        return False
    if flag2.coloured != {bijection[v] for v in flag1.coloured}:
        return False
    # Check edge mapping
    mapped_edges = set()
    for edge in flag1.edges:
        mapped_edge = frozenset(bijection[v] for v in edge)
        mapped_edges.add(mapped_edge)
    if mapped_edges != {frozenset(edge) for edge in flag2.edges}:
        return False
    return True


def are_isomorphic(flag1, flag2):
    """Checks if two flags are isomorphic.

    """
    # Check number of edges and coloured vertices
    if len(flag1.edges) != len(flag2.edges):
        return False
    if len(flag1.coloured) != len(flag2.coloured):
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
        new_flag = flag(G.size, G.edges, list(ftype), G.coloured)
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
    n, edges, coloured = typed_flag.size, typed_flag.edges, typed_flag.coloured
    type_size = len(typed_flag.ftype)
    numerator = 0
    for ftype in permutations(range(n), type_size):
        if are_isomorphic(typed_flag, flag(n, edges, list(ftype), coloured)):
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
                       coloured=[x[0] for x in (f[2][2][1] if
                                                len(f[2]) > 1 else [])])
                  for f in base_flags_0]

    factor_flags = []
    for a in factor_flags_0:
        b = []
        for c in a:
            b += [flag(size=c[0],
                       edges=[list(d) for d in c[2][0][1]],
                       ftype=[d for d in c[1]],
                       coloured=[x[0] for x in (c[2][2][1] if
                                                len(c[2]) > 1 else [])])]
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

    # Get list of valid 6-vertex flags, where the type is an edge, and the unlabeled
    # vertices induce an edge
    pre_processed = []
    for f in base_flags:
        # Keep only flags with at least two edges
        valid = len(f.edges) >= 2
        # Keep those with two disjoint edges
        valid = valid or next((True for a, b in combinations(f.edges, 2)
                               if len(set(a + b)) == 6), False)
        if valid:
            for ftype in permutations(range(6), 3):
                unlabeled_vertices = [v for v in range(6) if v not in ftype]
                e = [set(x) for x in f.edges]
                if set(ftype) in e and set(unlabeled_vertices) in e:
                    # If typed and untyped vertices induce edges
                    new_flag = flag(6, f.edges, list(ftype))
                    contained = next((True for other in pre_processed
                                      if are_isomorphic(new_flag, other)), False)
                    if not contained:
                        pre_processed += [new_flag]

    # Include only the special flags (in the same way as in the flag algebras
    # package)
    flag_terms = []
    for f in pre_processed:
        unlabeled_vertices = [v for v in range(6) if v not in f.ftype]
        A = 0
        for xx in unlabeled_vertices:
            touch = len(f.subflag(f.ftype + [xx]).edges) - 1
            if touch not in [0, 1]:
                A += 1
        if A == 2:
            flag_terms += [f]

    # Average the terms to produce a 0-quantum flag using the base flags
    objective = [0 for _ in range(len(base_flags))]
    for i, G in enumerate(base_flags):
        objective[i] = sum([average(term) for term in flag_terms
                            if are_isomorphic(G,
                                              flag(term.size, term.edges, coloured=term.coloured))])

    # Lower bound on average degree assumption (edge - 4641/10000 >= 0)
    avg_degree_minus_0_4641 = [
        (-Fraction(4641, 10000), flag(3, coloured=[0, 1, 2])),
        (-Fraction(4641, 10000), flag(3, coloured=[0, 1])),
        (-Fraction(4641, 10000), flag(3, coloured=[0])),
        (-Fraction(4641, 10000), flag(3, coloured=[])),
        (-Fraction(4641, 10000), flag(3, [[0, 1, 2]], coloured=[0, 1, 2])),
        (-Fraction(4641, 10000), flag(3, [[0, 1, 2]], coloured=[0, 1])),
        (-Fraction(4641, 10000), flag(3, [[0, 1, 2]], coloured=[0])),
        (1-Fraction(4641, 10000), flag(3, [[0, 1, 2]]))
    ]
    positivity_assumptions = [avg_degree_minus_0_4641]
    return objective, positivity_assumptions


def prop_3_3():
    """Returns the objective function as a quantum flag and the positivity
    assumptions as a list of quantum flags of Proposition 3.3.

    """
    objective = [
        (1, flag(3, edges=[[0, 1, 2]])),
        (1, flag(3, edges=[[0, 1, 2]], coloured=[1, 2])),
        (-(1 - Fraction(1, 4000)), flag(3, coloured=[2]))
    ]

    # Define assumptions

    # Local optimality
    Cp0_minus_Bp0 = [
        (1, flag(3, edges=[[0, 1, 2]], ftype=[0], coloured=[2])),
        (-1, flag(3, edges=[[0, 1, 2]], ftype=[0]))
    ]
    Cp1_minus_Bp1 = [
        (1, flag(3, edges=[[0, 1, 2]], ftype=[2], coloured=[2])),
        (-1, flag(3, edges=[[0, 1, 2]], ftype=[2], coloured=[1, 2]))
    ]

    # Degree regularity
    degeq_00 = [
        (1, flag(4, ftype=[2, 3], edges=[[0, 1, 2]], coloured=[0, 1])),
        (1, flag(4, ftype=[1, 3], edges=[[0, 1, 2]], coloured=[0])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[])),
        (1, flag(4, ftype=[2, 3], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[2, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[2, 1], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[])),
        (1, flag(4, ftype=[3, 2], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[3, 1], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[])),
        (-1, flag(4, ftype=[3, 2], edges=[[0, 1, 2]], coloured=[0, 1])),
        (-1, flag(4, ftype=[3, 1], edges=[[0, 1, 2]], coloured=[0])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[])),
        (-1, flag(4, ftype=[3, 2], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 2], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[1, 2], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[])),
        (-1, flag(4, ftype=[2, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[0, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[])),
    ]

    minus_degeq_00 = [(-coef, flag) for coef, flag in degeq_00]

    degeq_01 = [
        (1, flag(4, ftype=[3, 0], edges=[[1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[2, 0], edges=[[1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[2, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[2, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[2, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2]], coloured=[0, 1])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2]], coloured=[0])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[2, 0], edges=[[0, 1, 2], [0, 1, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0])),
        (-1, flag(4, ftype=[3, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[2, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0]))
    ]

    minus_degeq_01 = [(-coef, flag) for coef, flag in degeq_01]

    degeq_11 = [
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0, 1, 2, 3])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[1, 0], edges=[[1, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2, 3])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[2, 0], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (1, flag(4, ftype=[3, 0], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2, 3])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1, 2])),
        (1, flag(4, ftype=[0, 1], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[0, 1, 2, 3])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[0, 1], edges=[[1, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2, 3])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[0, 2], edges=[[0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 2, 3]], coloured=[0, 1])),
        (-1, flag(4, ftype=[0, 3], edges=[[0, 1, 3], [0, 2, 3], [1, 2, 3]], coloured=[0, 1, 2, 3])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1, 2])),
        (-1, flag(4, ftype=[1, 0], edges=[[0, 1, 2], [0, 1, 3], [0, 2, 3]], coloured=[0, 1]))
    ]

    minus_degeq_11 = [(-coef, flag) for coef, flag in degeq_11]

    # Best fit assumption (C - 2/5 >= 0)
    C_minus_2_5 = [
        (-Fraction(2, 5), flag(3, coloured=[0, 1, 2])),
        (-Fraction(2, 5), flag(3, coloured=[0, 1])),
        (-Fraction(2, 5), flag(3, coloured=[0])),
        (-Fraction(2, 5), flag(3, coloured=[])),
        (-Fraction(2, 5), flag(3, [[0, 1, 2]], coloured=[0, 1, 2])),
        (-Fraction(2, 5), flag(3, [[0, 1, 2]], coloured=[0, 1])),
        (Fraction(3, 5), flag(3, [[0, 1, 2]], coloured=[0])),
        (-Fraction(2, 5), flag(3, [[0, 1, 2]]))
    ]

    positivity_assumptions = [
        # Local optimality
        Cp0_minus_Bp0, Cp1_minus_Bp1,
        # Degree regularity
        degeq_00, minus_degeq_00,
        degeq_01, minus_degeq_01,
        degeq_11, minus_degeq_11,
        # Best fit assumption
        C_minus_2_5
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
  $ python validator.py prop_3_1_c5k4.pickle
""")
parser.add_argument("certificate", type=str, help='path of the pickle file',
                    choices=['prop_3_1_c5k4.pickle', 'prop_3_1_c7.pickle',
                             'prop_3_2_c5k4.pickle', 'prop_3_2_c7.pickle',
                             'prop_3_3_c5k4.pickle', 'prop_3_3_c7.pickle'])

certificate_path = parser.parse_args().certificate

(base_flags,
 target,
 factor_flags,
 positives) = recover_from_certificate(certificate_path)

# Obtain objective and assumptions
# =================================================
if certificate_path in ['prop_3_1_c7.pickle', 'prop_3_1_c5k4.pickle']:
    objective, positivity_assumptions = prop_3_1()
elif certificate_path in ['prop_3_2_c5k4.pickle', 'prop_3_2_c7.pickle']:
    objective, positivity_assumptions = prop_3_2(base_flags)
elif certificate_path in ['prop_3_3_c5k4.pickle', 'prop_3_3_c7.pickle']:
    objective, positivity_assumptions = prop_3_3()

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
    transposed = list(map(list, zip(*densities)))
    n_prev_blocks = sum([len(factor_flags[idx]) for idx in range(idx)])
    start, end = n_prev_blocks, n_prev_blocks + len(factor_flags[idx])
    if transposed == positives[start:end]:
        print("Positivity assumption %d matches." % (idx + 1))
    else:
        print("ERROR: POSITIVITY ASSUMPTION %d DOES NOT MATCH!" % (idx + 1))
        exit(1)

print("\nCOMPLETED: THE OBJECTIVES AND ASSUMPTIONS CORRESPOND TO THOSE IN THE PROPOSITION!")
