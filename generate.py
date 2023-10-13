import re
from sympy import Matrix, zeros
from sympy.polys.matrices import DomainMatrix
from itertools import groupby

import reference

HYBRIDIZATION = False

species_to_solve = ['C18H14', (0,0,16,2)]

re_elems = re.compile(r'([A-Z]+?)(?=\d+)?')
re_with_counts = re.compile(r'([A-Z]+?)(\d+)?')

elem_num = {el:i for i, el in enumerate(set([e for r in reference.reference.keys() for e in re_elems.findall(r.split('_')[-2]) if not (HYBRIDIZATION and 'C' in e)]))}
if HYBRIDIZATION: 
	l = len(elem_num)
	elem_num.update({f'C{i}':l+i for i in range(4)})

dim = len(elem_num)
num_elem = dict(map(reversed, elem_num.items()))

def str_to_vec(s, hyb=None):
	r_vec = zeros(dim, 1)
	for el, n in re_with_counts.findall(s):
		if HYBRIDIZATION and el=='C':
			for i, c in enumerate(hyb):
				r_vec[elem_num[f'C{i}']] = int(c)
		else:
			r_vec[elem_num[el]] = int(n) if n else 1
	return r_vec

basis = [str_to_vec(*species_to_solve)] + [str_to_vec(r.split('_')[-2], reference.hybridization_1[r] if HYBRIDIZATION else None) for r in reference.reference.keys()]
basis_l = len(basis)

def minimal_nullspace(A):
	n = A.to_field().nullspace().to_Matrix()
	found = n.shape[0]>=1

	return (n if n.shape[0]==1 and not 0 in n else None), found

S = []
vec_nums = [0, 1]
i = 2
while len(vec_nums) > 1:
	A = DomainMatrix.from_Matrix(Matrix.hstack(*(basis[k] for k in vec_nums)))
	n, found = minimal_nullspace(A)

	if found:
		if n is not None:
			S.append([vec_nums.copy(), n/n[0]])
		vec_nums.pop()

	if i < basis_l: 
		vec_nums.append(i)
	else:
		while vec_nums[-1]>=basis_l and len(vec_nums)>1:
			vec_nums.pop()
		i = vec_nums.pop() + 1
		if i < basis_l:
			vec_nums.append(i)

	i += 1

def coef_to_str(c):
	return str(c) if c!=1 else ''

def vector_to_str(v):
	s = ''

	if HYBRIDIZATION:
		c_sum = sum(v[n] for (e, n) in elem_num.items() if 'C' in e)
		if c_sum > 0:
			s += ''.join(f'C{coef_to_str(c_sum)}')

	for i in range(dim):
		if (HYBRIDIZATION and 'C' in num_elem[i]): continue 
		if v[i] > 0: s += ''.join(f'{num_elem[i]}{coef_to_str(v[i])}')
	return s

def reaction_to_str(rxn, cfs):
	reagents = [[r, c] for (r, c) in zip(rxn, cfs) if c>0]
	products = [[r, c] for (r, c) in zip(rxn, cfs) if c<0]

	s = ''
	s += '+'.join('*'.join(filter(None, [coef_to_str(c), vector_to_str(basis[r])])) for (r, c) in reagents)
	s += '->'
	s += '+'.join('*'.join(filter(None, [coef_to_str(-c), vector_to_str(basis[r])])) for (r, c) in products)

	return s

def is_integer_vec(v):
	for i in range(1,len(v)):
		if not v[i].is_integer:
			return False
	return True

def sort_func(rxn):
	return sum(abs(rxn[1]))+abs(sum(rxn[1]))

n_rxn = len(S)
S = list(filter(lambda rxn: is_integer_vec(rxn[1]), S))
S = sorted(S, key=sort_func)

layers = [list(g[1]) for g in groupby(S, sort_func)]

print(f'Reactions found: {n_rxn}, integer: {len(S)}')
for i in range(len(layers)):
	print(f'Layer {i}:')
	print(*(reaction_to_str(*rxn) for rxn in layers[i]), sep='\n')
