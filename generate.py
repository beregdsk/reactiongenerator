import re
from sympy import QQ
from sympy.polys.matrices import DomainMatrix
from itertools import groupby

import reference

class ReactionGenerator:
    def __init__(self, species_to_solve, ref, hyb=False):
        self.species_to_solve = species_to_solve
        self.hybridization = hyb
        
        self.re_elems = re.compile(r'([A-Z]+?)(?=\d+)?')
        self.re_with_counts = re.compile(r'([A-Z]+?)(\d+)?')
        elements = sorted(set([e for r in ref.keys() 
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
        r_vec = [0]*self.dim
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

    def generate(self, fraction=None):
        if fraction is not None:
            start = int(fraction[0]/fraction[1]*(len(self.basis)-1) + 1)
            end = int((fraction[0]+1)/fraction[1]*(len(self.basis)-1) + 1)
        else:
            start = 1
            end = len(self.basis)
        
        S = []
        vec_nums = [0, start]
        while len(vec_nums) > 1 and vec_nums[1] < end:
            A = DomainMatrix.from_list([self.basis[k] for k in vec_nums], QQ).transpose()
            n, found = self.minimal_nullspace(A)

            if found:
                if n is not None:
                    S.append([vec_nums.copy(), n/n[0]])

                last = vec_nums.pop()
                if last+1 >= len(self.basis):
                    last = vec_nums.pop()
                vec_nums.append(last+1)
            else:
                if vec_nums[-1]+1 < len(self.basis):
                    vec_nums.append(vec_nums[-1]+1)
                else:
                    last = vec_nums.pop()
                    i = 1
                    while last == len(self.basis)-i:
                        if len(vec_nums) <= 1: 
                            return S
                        last = vec_nums.pop()
                        i += 1
                    vec_nums.append(last+1)

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
    from multiprocessing import Pool

    MP = False

    gen = ReactionGenerator(('C18H14', (0,16,2)), reference.hybridization, False)
    if MP:
        rxns = []

        n_proc = 4
        with Pool(n_proc) as p:
            for res in p.imap_unordered(gen.generate, ((i, n_proc) for i in range(n_proc))):
                rxns.extend(res)
    else:
        rxns = gen.generate()

    int_rxns = list(filter(lambda rxn: gen.is_integer_vec(rxn[1]), rxns))
    int_rxns = sorted(int_rxns, key=sort_func)

    layers = [list(g[1]) for g in groupby(int_rxns, sort_func)]

    print(f'Reactions found: {len(rxns)}, integer: {len(int_rxns)}')
    input()
    for i in range(len(layers)):
        print(f'Layer {i+1}:')
        print(*(gen.reaction_to_str(*rxn) for rxn in layers[i]), sep='\n')
