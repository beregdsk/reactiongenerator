import re
from sympy import Matrix, zeros
from sympy.polys.matrices import DomainMatrix
from itertools import groupby
from collections import deque

import reference

class ReactionGenerator:
    def __init__(self, species_to_solve, ref, hyb=False):
        self.species_to_solve = species_to_solve
        self.hybridization = hyb
        
        self.re_elems = re.compile(r'([A-Z]+?)(?=\d+)?')
        self.re_with_counts = re.compile(r'([A-Z]+?)(\d+)?')
        elements = sorted(set([e for r in reference.reference.keys() 
                          for e in self.re_elems.findall(r.split('_')[-2]) 
                          if not (self.hybridization and 'C' in e)]))
        self.elem_num = {el:i for i, el in enumerate(elements)}
        
        if hyb: self.elem_num.update({f'C{i}':len(self.elem_num)+i for i in range(4)})

        self.num_elem = dict(map(reversed, self.elem_num.items()))

        self.dim = len(self.elem_num)
        
        self.basis = [self.str_to_vec(*species_to_solve)] + \
            [self.str_to_vec(r.split('_')[-2], 
                ref[r] if self.hybridization else None) 
                for r in ref.keys()]

    def str_to_vec(self, s, hyb=None):
        r_vec = zeros(self.dim, 1)
        for el, n in self.re_with_counts.findall(s):
            if self.hybridization and el=='C':
                for i, c in enumerate(hyb):
                    r_vec[self.elem_num[f'C{i}']] = int(c)
            else:
                r_vec[self.elem_num[el]] = int(n) if n else 1
        return r_vec

    @staticmethod
    def minimal_nullspace(A):
        n = A.to_field().nullspace().to_Matrix()
        found = n.shape[0]>=1

        return (n if n.shape[0]==1 and not 0 in n else None), found

    def generate(self):
        S = []
        vec_nums = deque([0, 1])
        i = len(vec_nums)
        while len(vec_nums) > 1:
            A = DomainMatrix.from_Matrix(Matrix.hstack(*([self.basis[k] for k in vec_nums])))
            n, found = self.minimal_nullspace(A)

            if found:
                if n is not None:
                    S.append([vec_nums.copy(), n/n[0]])
                vec_nums.pop()

            if i < len(self.basis): 
                vec_nums.append(i)
            else:
                while vec_nums[-1]>=len(self.basis) and len(vec_nums)>1:
                    vec_nums.pop()
                i = vec_nums.pop() + 1
                if i < len(self.basis):
                    vec_nums.append(i)

            i += 1

        return S

    @staticmethod
    def coef_to_str(c):
        return str(c) if c!=1 else ''

    def vector_to_str(self, v):
        s = ''

        if self.hybridization:
            c_sum = sum(v[n] for (e, n) in self.elem_num.items() if 'C' in e)
            if c_sum > 0:
                s += ''.join(f'C{self.coef_to_str(c_sum)}')

        for i in range(self.dim):
            if (self.hybridization and 'C' in self.num_elem[i]): continue 
            if v[i] > 0: s += ''.join(f'{self.num_elem[i]}{self.coef_to_str(v[i])}')
        return s

    def reaction_to_str(self, rxn, cfs):
        reagents = [[r, c] for (r, c) in zip(rxn, cfs) if c>0]
        products = [[r, c] for (r, c) in zip(rxn, cfs) if c<0]

        s = ''
        s += '+'.join(
                '*'.join(
                    filter(None, [self.coef_to_str(c), self.vector_to_str(self.basis[r])])) 
                for (r, c) in reagents)
        s += '->'
        s += '+'.join(
                '*'.join(
                    filter(None, [self.coef_to_str(-c), self.vector_to_str(self.basis[r])])) 
                for (r, c) in products)

        return s

    @staticmethod
    def is_integer_vec(v):
        for i in range(1,len(v)):
            if not v[i].is_integer:
                return False
        return True

def sort_func(rxn):
    return sum(abs(rxn[1]))+abs(sum(rxn[1]))

if __name__ == '__main__':
    gen = ReactionGenerator(('C18H14', (0,0,16,2)), reference.hybridization)

    rxns = gen.generate()
    int_rxns = list(filter(lambda rxn: gen.is_integer_vec(rxn[1]), rxns))
    int_rxns = sorted(int_rxns, key=sort_func)

    layers = [list(g[1]) for g in groupby(int_rxns, sort_func)]

    print(f'Reactions found: {len(rxns)}, integer: {len(int_rxns)}')
    for i in range(len(layers)):
        print(f'Layer {i+1}:')
        print(*(gen.reaction_to_str(*rxn) for rxn in layers[i]), sep='\n')
