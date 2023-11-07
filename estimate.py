from multiprocessing import Process, Manager
from functools import partial
import numpy as np
import ctypes
import sys

import reference

rgen = ctypes.CDLL('./rgen{}'.format('.dll' if sys.platform[:3].lower() == 'win' else '.so'))
rgen.generate.argtypes = (np.ctypeslib.ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
                          ctypes.POINTER(ctypes.c_int),
                          np.ctypeslib.ndpointer(ctypes.c_int, flags='C_CONTIGUOUS'),
                          ctypes.c_int, ctypes.c_int, ctypes.c_int)
rgen.generate.restype = None

path_to_files = 'PORPHYRIN_FINAL_NOSO/'

ehtocal = 627.51

class EnthalpyEstimator:
    def __init__(self, filename, reference, isodesmic=True):
        self.filename = filename
        self.reference = reference
        self.isodesmic = isodesmic

        self.init_ref()
        self.init_gen()

    def init_ref(self):
        self.ref = {}
        
        def read_energies(entry, filename):
            try:
                with open(path_to_files + filename + 'PBE0_ccsdt_cc-pwcvqz_ecp.txt', 'r') as f:
                    for l in f.readlines()[2:]:
                        if '=' in l:
                            l_split = [s.strip() for s in l.split('=')]
                            entry.setdefault('E_all', {})[l_split[0]] = float(l_split[-1])
            except Exception:
                print('Error reading file ' + path_to_files + filename + 'wB97XD_ccsdt_cc-pvqz_ecp.txt')

        def read_Hcorr(entry, filename):
            try:
                with open(path_to_files + filename + 'PBE0.td', 'r') as f:
                    entry.setdefault('E_all', {})['Hcorr'] = [float(l.split()[-1]) for l in f
                                    if 'Thermal correction to Enthalpy' in l][0]
            except Exception:
                print('Error reading file ' + path_to_files + filename + 'wB97XD.td')

        def read_atoms(entry, filename):
            try:
                with open(path_to_files + filename + 'PBE0.xyz', 'r') as f:
                    atoms = {}
                    for l in f:
                        if ' ' not in l: continue
                        a = l.split()[0].strip()
                        atoms[a] = atoms.get(a, 0) + 1

                    entry['atoms'] = atoms
            except Exception:
                print('Error reading file ' + path_to_files + filename + 'PBE0.xyz')

        for r in self.reference:
            entry = {}
            
            read_energies(entry, r)
            read_Hcorr(entry, r)
            read_atoms(entry, r)

            entry['E_all']['dfH'] = self.reference[r][0]
            entry['E_all']['u'] = self.reference[r][1]

            self.ref[r] = entry

        entry = {}
        
        read_energies(entry, self.filename)
        read_Hcorr(entry, self.filename)
        read_atoms(entry, self.filename)

        self.ref[self.filename] = entry

    def init_gen(self):
        self.elements = np.array(sorted({a for r in self.ref.values() for a in r['atoms'].keys()}))

        self.dim = len(self.elements)
        
        ref_ordered = [self.ref[self.filename]] + [self.ref[r] for r in self.ref if r != self.filename]
        self.basis = np.array([[r['atoms'].get(e, 0) for e in self.elements] for r in ref_ordered], dtype=np.int32)

    def run(self, max_rxns=1000):
        S = np.empty((max_rxns, len(self.basis)), dtype=np.int32)
        len_S = ctypes.pointer(ctypes.c_int(0))
        
        rgen.generate(S, len_S, self.basis, self.dim, len(self.basis), max_rxns)
        rxns = S[:int(len_S.contents.value)]
        # if len(rxns) > 0:
        #     rxns = rxns[np.apply_along_axis(lambda row: abs(np.sum(row)) + np.sum(np.abs(row)), 1, rxns).argsort()]

        self.reactions = rxns
        
        # filenames = [self.filename] + [name for name in self.ref 
        #                              if name != self.filename]

        # with open(path_to_files + self.filename + 'wB97XD.all', 'w') as f:
        #     for rxn in self.reactions:
        #         entry = []
        #         entry.append(self.reaction_to_str(rxn))
        #         entry.append(sum(rxn))
        #         entry.append(sum(abs(rxn)))
        #         entry.append(abs(entry[1]+entry[2]))

        #         E_sum = {}
        #         reactants = np.where(rxn != 0)[0]
        #         for r in reactants:
        #             c_frac = rxn[r].item()/rxn[reactants[0]].item()
        #             E_all = self.ref[filenames[r]]['E_all']

        #             for n, e in E_all.items():
        #                 if n=='u':
        #                     E_sum[n] = E_sum.get(n, 0) + (c_frac*e)**2
        #                 else:
        #                     E_sum[n] = E_sum.get(n, 0) + c_frac*e
                
        #         E_CV = E_sum['E(CV correction (CBS))']
        #         E_IT = E_sum['E(DLPNO-CCSD(T)/TZ-IT)']-E_sum['E(DLPNO-CCSD(T)/CC-PVTZ)']
        #         E_SO = E_sum['E(SO_EFFECT_TZ)']
        #         E_TPNO = E_sum['E(DLPNO_TPNO/CC-PVTZ))']-E_sum['E(DLPNO-CCSD(T)/CC-PVTZ)']

        #         E_r = E_sum['E(DLPNO-CCSD(T)/CBS)'] + E_SO + E_IT + E_TPNO + E_CV
        #         E_r_HF = E_sum['E(HF/CBS)'] + E_SO
        #         E_corr_r = E_r - E_r_HF
                
        #         H_r = E_r + E_sum['Hcorr']
        #         H_f = E_sum['dfH'] - H_r*ehtocal

        #         u = E_sum['u']**0.5

        #         entry.append(round(H_f, 3))
        #         entry.append(round(E_corr_r*ehtocal, 3))
        #         entry.append(round(E_sum['Hcorr']*ehtocal, 3))
        #         entry.append(entry[-1]+entry[-2])
        #         entry.append(round(u, 3))
        #         entry.append(round(2*u, 3))
        #         entry.append(round(H_r*ehtocal, 3))
        #         entry.append(0)
        #         entry.append(round(E_CV*ehtocal, 3))
        #         entry.append(round(E_IT*ehtocal, 3))
        #         entry.append(round(E_TPNO*ehtocal, 3))
        #         entry.append(round(E_SO*ehtocal, 3))

        #         entry.append(self.reaction_to_str(rxn))

        #         f.write(';'.join(map(str, entry)) + '\n')

    @staticmethod
    def coef_to_str(c):
        return str(c) if c!=1 else ''

    def vector_to_str(self, v):
        els = np.where(v != 0)[0]
        el_n = np.column_stack([self.elements[els], v[els]])
        el_n[el_n[:, 1] == '1', 1] = ''
        s = ''.join(''.join(el_n.flatten()))
                    
        return s

    def reaction_to_str(self, rxn):
        r = rxn * np.sign(rxn[0])
        reagents = np.where(r > 0)[0]
        products = np.where(r < 0)[0]

        s = ''
        s += '+'.join(
                ''.join([self.coef_to_str(r[i].item()), self.vector_to_str(self.basis[i])]) 
                for i in reagents)
        s += '->'
        s += '+'.join(
                ''.join([self.coef_to_str(-r[i].item()), self.vector_to_str(self.basis[i])]) 
                for i in products)

        return s

if __name__ == '__main__':
    import timeit
    from tqdm import tqdm
    import matplotlib.pyplot as plt

    times = []
    for i in tqdm(range(2, len(reference.reference))):
        subdict = dict(list(reference.reference.items())[:i])
        est = EnthalpyEstimator('m2_00_Ag_model6_D4h_plane-', subdict, isodesmic=False)

        times.append(timeit.Timer(lambda: est.run(max_rxns=50000)).timeit(number=3)/3)

    plt.plot(range(2, len(times)+2), times)
    plt.show()