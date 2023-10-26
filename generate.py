from math import lcm
from sympy import QQ, sign
from sympy.polys.matrices import DomainMatrix
from itertools import groupby
from rdkit import Chem
from rdkit.Chem import rdqueries

import reference

class ReactionGenerator:
    def __init__(self, species_to_solve, smi, isodesmic=True):
        self.isodesmic = isodesmic
        self.smiles = [species_to_solve] + [s for s in smi if s != '']
        self.mols = {s:Chem.AddHs(Chem.MolFromSmiles(s)) for s in self.smiles}
        
        atomic_nums = sorted(set([atom.GetAtomicNum() for mol in self.mols.values() 
                                for atom in mol.GetAtoms()]), reverse=True)
        self.elem_num = {el:i for i, el in enumerate(atomic_nums)}
        self.bonds = [Chem.MolFromSmarts(s) for s in ('[C]-[H]', '[C]-[C]', '[C]=[C]', '[C]#[C]', '[C]:[C]')]

        self.dim = len(self.elem_num) + (len(self.bonds) if isodesmic else 0)
        
        self.basis = [self.smi_to_vec(s) for s in self.smiles]

    def smi_to_vec(self, s):
        r_vec = [0]*self.dim
        
        for el, i in self.elem_num.items():
            q = rdqueries.AtomNumEqualsQueryAtom(el)
            r_vec[i] = len(self.mols[s].GetAtomsMatchingQuery(q))

        if self.isodesmic:
            for i, b in enumerate(self.bonds):
                r_vec[len(self.elem_num)+i] = len(self.mols[s].GetSubstructMatches(b))

        return r_vec

    @staticmethod
    def minimal_nullspace(A):
        n = A.to_field().nullspace().to_Matrix()
        found = n.shape[0]>=1
        
        return (n if n.shape[0]==1 and not 0 in n else None), found

    def generate(self, fraction=(0,1), max_found=1000):
        start = int(fraction[0]/fraction[1]*(len(self.basis)-1) + 1)
        end = int((fraction[0]+1)/fraction[1]*(len(self.basis)-1) + 1)
        
        S = []
        vec_nums = [0, start]
        while len(vec_nums) > 1 and vec_nums[1] < end and len(S) < max_found/fraction[1]:
            A = DomainMatrix.from_list([self.basis[k] for k in vec_nums], QQ).transpose()
            n, found = self.minimal_nullspace(A)

            if found:
                if n is not None:
                    l = lcm(*(a.q for a in n))
                    S.append([vec_nums.copy(), n*l*sign(n[0])])

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

        pt = Chem.GetPeriodicTable()
        for el, i in self.elem_num.items():
            if v[i] > 0: s += ''.join(f'{pt.GetElementSymbol(el)}{self.coef_to_str(v[i])}')
        return s

    def reaction_to_str(self, rxn, cfs):
        reagents = [[r, c] for (r, c) in zip(rxn, cfs) if c>0]
        products = [[r, c] for (r, c) in zip(rxn, cfs) if c<0]

        s = ''
        s += '+'.join(
                ''.join(
                    filter(None, [self.coef_to_str(c), self.vector_to_str(self.basis[r])])) 
                for (r, c) in reagents)
        s += '->'
        s += '+'.join(
                ''.join(
                    filter(None, [self.coef_to_str(-c), self.vector_to_str(self.basis[r])])) 
                for (r, c) in products)

        return s
    
    def reaction_to_smi(self, rxn, cfs):
        reagents = [[r, c] for (r, c) in zip(rxn, cfs) if c>0]
        products = [[r, c] for (r, c) in zip(rxn, cfs) if c<0]

        s = ''
        s += '+'.join(
                ''.join(
                    filter(None, [self.coef_to_str(c), self.smiles[r]])) 
                for (r, c) in reagents)
        s += '->'
        s += '+'.join(
                ''.join(
                    filter(None, [self.coef_to_str(-c), self.smiles[r]])) 
                for (r, c) in products)

        return s
    
    def is_isodesmic(self, rxn, cfs):
        smarts = ['[C]-[H]', '[C]-[C]', '[C]=[C]', '[C]#[C]', '[C]:[C]']
        smarts_mols = [Chem.MolFromSmarts(s) for s in smarts]
        t_s = [0]*len(smarts)
        for r, c in zip(rxn, cfs):
            mol = self.mols[self.smiles[r]]
            for i, m in enumerate(smarts_mols):
                t_s[i] += c*len(mol.GetSubstructMatches(m))

        s = sum(abs(x) for x in t_s)
        return s == 0
    
    def is_homodesmotic(self, rxn, cfs):
        smarts = ['[^3]-[^3]', '[^2]-[^3]', '[^1]-[^3]', '[^2]-[^2]', '[^2]-[^1]', '[^1]-[^1]', '[^2]=[^2]', '[^2]:[^2]', '[^1]#[^1]',
                '[C^3H3]', '[C^3H2]', '[C^3H]', '[C^3H0]', '[C^2H2]', '[C^2H]', '[C^2H0]', '[C^1H]', '[C^1H0]']
        smarts_mols = [Chem.MolFromSmarts(s) for s in smarts]
        t_s = [0]*len(smarts)
        for r, c in zip(rxn, cfs):
            mol = self.mols[self.smiles[r]]
            for i, m in enumerate(smarts_mols):
                t_s[i] += c*len(mol.GetSubstructMatches(m))

        s = sum(abs(x) for x in t_s)
        return s == 0
    
    def is_hypohomodesmotic(self, rxn, cfs):
        smarts = ['[^3]', '[^2]', '[^1]',
                '[C^3H3]', '[C^3H2]', '[C^3H]', '[C^3H0]', '[C^2H2]', '[C^2H]', '[C^2H0]', '[C^1H]', '[C^1H0]']
        smarts_mols = [Chem.MolFromSmarts(s) for s in smarts]
        t_s = [0]*len(smarts)
        for r, c in zip(rxn, cfs):
            mol = gen.mols[self.smiles[r]]
            for i, m in enumerate(smarts_mols):
                t_s[i] += c*len(mol.GetSubstructMatches(m))

        s = sum(abs(x) for x in t_s)
        return s == 0


def sort_func(rxn):
    return sum(abs(rxn[1]))+abs(sum(rxn[1]))

if __name__ == '__main__':
    from multiprocessing import Pool

    MP = True

    gen = ReactionGenerator('c1cc2cccc3c2c4c1cccc4cc3', reference.smiles, True)

    if MP:
        rxns = []

        n_proc = 4
        with Pool(n_proc) as p:
            for res in p.imap_unordered(gen.generate, ((i, n_proc) for i in range(n_proc))):
                rxns.extend(res)
    else:
        rxns = gen.generate()

    rxns = sorted(rxns, key=sort_func)
    layers = [list(g[1]) for g in groupby(rxns, sort_func)]

    print(f'Reactions found: {len(rxns)}')
    input()
    for i in range(int(len(layers)/10)):
        print(f'Layer {i+1}:')
        print(*(gen.reaction_to_str(*rxn) for rxn in layers[i]), sep='\n')