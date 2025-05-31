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
# Global variables:
# verbose
# typed_flags
# type_matrices
# positives

def is_positive_semidefinite(M):
    """Performs an LDL decomposition of the matrix to verify it is positive
    semidefinite.

    """
    n = len(M)
    # Check if the matrix is square
    for row in M:
        if len(row) != n:
            return False
    # Check symmetry
    for i in range(n):
        for j in range(i + 1, n):
            if M[i][j] != M[j][i]:
                return False
    # Initialize L and D
    L = [[Fraction(0) for _ in range(n)] for __ in range(n)]
    D = [Fraction(0) for _ in range(n)]
    for i in range(n):
        L[i][i] = Fraction(1)
    # Compute L and D
    for j in range(n):
        sum_prev = sum(L[j][k] ** 2 * D[k] for k in range(j))
        D[j] = M[j][j] - sum_prev
        if D[j] == 0:
            for i in range(j + 1, n):
                L[i][j] = Fraction(0)
        else:
            for i in range(j + 1, n):
                sum_prev_i = sum(L[i][k] * L[j][k] * D[k] for k in range(j))
                L[i][j] = (M[i][j] - sum_prev_i) / D[j]
    # Check D positivity
    if any(d < 0 for d in D):
        return False
    # Verify M = LDL^T
    for i in range(n):
        for k in range(n):
            total = sum(L[i][j] * D[j] * L[k][j] for j in range(n))
            if total != M[i][k]:
                return False
    return True


def inner_product(M1, M2):
    """Calculates the inner product of two matrices.

    """
    sum = 0
    for i in range(len(M1)):
        for j in range(len(M1)):
            sum += M1[i][j] * M2[i][j]
    return sum


def reconstruct_matrix(size, flat):
    """Given a flattened matrix, returns the matrix as a list of rows.

    """
    matrix = [[0] * size for _ in range(size)]
    for row in range(size):
        start_idx = row * (2 * size - row + 1) // 2
        for k in range(row, size):
            value = flat[start_idx + (k - row)]
            matrix[row][k] = matrix[k][row] = value
    return matrix


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


def get_multiplication_matrix(flags, G):
    """Given a set of flags from the same type, calculates the multiplication matrix
    for a 0-flag G.

    """
    n = G.size
    type_size = len(flags[0].ftype)
    flag_size = flags[0].size
    M = [[Fraction(0) for _ in range(len(flags))] for _ in range(len(flags))]
    V = range(n)
    for ftype in permutations(V, type_size):
        new_flag = flag(G.size, G.edges, list(ftype), G.coloured)
        unlabeled_vertices = [v for v in V if v not in ftype]
        for V1 in combinations(unlabeled_vertices, flag_size - type_size):
            flag1 = new_flag.subflag(ftype + V1)
            for idx, f in enumerate(flags):
                if are_isomorphic(flag1, f):
                    remaining_vertices = [v for v in unlabeled_vertices if v not in V1]
                    for V2 in combinations(remaining_vertices, flag_size - type_size):
                        flag2 = new_flag.subflag(ftype + V2)
                        for idx2, g in enumerate(flags):
                            if are_isomorphic(flag2, g):
                                M[idx][idx2] += 1
    denominator = comb(n, type_size) * factorial(type_size)
    denominator *= comb(n - type_size, flag_size - type_size)
    denominator *= comb(n - flag_size, flag_size - type_size)
    for i in range(len(flags)):
        for j in range(len(flags)):
            M[i][j] /= denominator
    return M


def get_inner_sum(F):
    """Given a 0-flag F, calculates the sum of inner products Σ<M_F,X>.

    """
    if verbose >= 2:
        print("\nSum of inner products:")
    inner_sum = 0
    for type_index, same_type_flags in enumerate(typed_flags):
        if verbose >= 2:
            print("Considering type", type_index)
        M_type = get_multiplication_matrix(same_type_flags, F)
        X_type = type_matrices[type_index]
        product = inner_product(M_type, X_type)
        inner_sum += product
        if verbose >= 3:
            print("\nMultiplication matrix for type %d:" % type_index)
            for row in M_type:
                print("[", *row, "]")
            print("X matrix for type %d:" % type_index)
            for row in X_type:
                print("[", *row, "]")
        if verbose >= 2:
            print("<M,X> =", product)
    return inner_sum


def get_positivity_sum(base_graph_index):
    """Given the index of a 0-flag F, calculates the positivity sum Σe*d(f, F).

    """
    if verbose >= 2:
        print("\nPositivity sum:")
    positivity_sum = 0
    for positivity_index, positivity_coefficients in enumerate(positives):
        e_coefficient = e_vector[positivity_index]
        pos_density_in_F = positivity_coefficients[base_graph_index]
        product = e_coefficient * pos_density_in_F
        positivity_sum += product
        if verbose >= 2:
            print("(", e_coefficient, ") (", pos_density_in_F, ") =", product)
    return positivity_sum


def recover_from_certificate(file_name):
    """Given a certificate, recovers the stored objects.

    """
    with open(file_name, 'rb') as handle:
        certificate = pickle.load(handle)
        base_flags_0 = certificate['base flags']
        typed_flags_0 = certificate['typed flags']
        types = list(typed_flags_0.keys())
        x_matrices = certificate['X matrices']
        target = certificate['target']
        positives = certificate['positives']
        e_vector = certificate['e vector']
        result = certificate['result']
        slack_vector = certificate['slack vector']
        maximize = certificate['maximize']

    # # Minimization changes the sign
    # if not maximize:
    #     result = -result

    # Get list of base flags and typed flags in format
    base_flags = [flag(size=f[0],
                       edges=[list(e) for e in f[2][0][1]], ftype=[],
                       coloured=[x[0] for x in (f[2][2][1] if
                                                len(f[2]) > 1 else [])])
                  for f in base_flags_0]
    typed_flags = []
    for ftype in types:
        same_type_flags = [flag(tflag[0],
                                [list(e) for e in tflag[2][0][1]],
                                list(tflag[1]),
                                [x[0] for x in (tflag[2][2][1]
                                                if len(tflag[2]) > 1 else [])])
                           for tflag in typed_flags_0[ftype]]
        typed_flags += [same_type_flags]
    return base_flags, typed_flags, x_matrices, target, positives, e_vector, result, slack_vector, maximize


# ---------------------------- Main ----------------------------

# Read certificate and ask for input
# =================================================
parser = argparse.ArgumentParser(
    description="A simple flag algebras certificate verifier.",
    formatter_class=argparse.RawTextHelpFormatter,
    epilog="""
examples:
  $ python verify.py c5k4-bad-missing.pickle

  $ python verify.py c5k4-bad-missing.pickle -p
  Enter two integers in the range [1,28080] (e.g., `1 14040'): 1 10

  $ python verifier.py c5k4-bad-missing.pickle -p -v 2
  Enter two integers in the range [1,28080] (e.g., `1 14040'): 5 5
""")
parser.add_argument("certificate", type=str, help='path of the pickle file')
parser.add_argument("-p", "--partial", action='store_true',
                    help='partial verification\n')
parser.add_argument("-v", "--verbose", type=int, choices=[1, 2, 3], default=1,
                    help="verbosity level\ndefault: 1")
args = parser.parse_args()

certificate_path = args.certificate
verbose = args.verbose
partial_verification = args.partial

(base_flags,
 typed_flags,
 x_matrices,
 target,
 positives,
 e_vector,
 result,
 slack_vector,
 maximize) = recover_from_certificate(certificate_path)

start, end = 1, len(base_flags)
if partial_verification:
    start, end = [int(x) for x in input("Enter two integers in the range [1,%d] (e.g., `1 %d'): " % (end, (end + 1) // 2)).split(" ")]


# Positive semidefinite verification
# =================================================
print("\nVerifying positive semidefiniteness")
type_matrices = []
for index, flat_matrix in enumerate(x_matrices):
    matrix = reconstruct_matrix(len(typed_flags[index]), flat_matrix)
    type_matrices += [matrix]
    if is_positive_semidefinite(matrix):
        print("Matrix %d is positive semidefinite" % (index + 1))
    else:
        print("VERIFICATION FAILED: MATRIX %d IS NOT POSITIVE SEMIDEFINITE!" % (index + 1))
        exit(1)
del x_matrices


# Constraint verification
# =================================================
print("\nVerifying constraints")
for index in range(start - 1, end):
    if verbose >= 2:
        print("--------------------------------------------------------------")
    print("Verifying flag %d" % (index + 1))

    F = base_flags[index]
    obj_fun_density = target[index]
    inner_sum = get_inner_sum(F)
    positivity_sum = get_positivity_sum(index)
    slack = slack_vector[index]

    if verbose >= 2:
        print("Fixed base graph: size: %d edges: %s ftype: %s coloured: %s" %
              (F.size, str(F.edges), str(F.ftype), str(list(F.coloured))))
        print("\nObjective function density:", obj_fun_density)
        print("Inner sum:", inner_sum)
        print("Positivity sum:", positivity_sum, '\n')

    if maximize:
        print("bound - d(T, F_%d) = Σ<M,X> + Σe*d(f, F_%d) + slack" % (index + 1, index + 1))
        print(result, '-', obj_fun_density, '=', inner_sum, '+', positivity_sum, '+', slack)
    else:
        print("d(T, F_%d) - bound = Σ<M,X> + Σe*d(f, F_%d) + slack" % (index + 1, index + 1))
        print(obj_fun_density, '-', result, '=', inner_sum, '+', positivity_sum, '+', slack)

    if slack < 0:
        print("VERIFICATION FAILED: SLACK IS NEGATIVE!")
        exit(1)

    sign = 1 if maximize else -1
    if sign*(result - obj_fun_density) == inner_sum + positivity_sum + slack:
        print("Equality verified\n")
    else:
        print("VERIFICATION FAILED: EQUALITY DOES NOT HOLD!")
        exit(1)

print("VERIFICATION COMPLETED SUCCESSFULLY!")
