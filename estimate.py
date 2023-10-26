from multiprocessing import Pool

from generate import ReactionGenerator
import reference

path_to_files = 'PAH_C32_ALL/'

ehtocal = 627.51

class EnthalpyEstimator:
    def __init__(self, filename, smi, name=None, isodesmic=True):
        self.filename = filename
        self.smi = smi
        self.name = name
        self.isodesmic = isodesmic

        self.init_ref()
        self.init_gen()

    def init_ref(self):
        self.ref = {}
        
        def read_energies(entry, filename):
            try:
                with open(path_to_files + filename + 'wB97XD_ccsdt_cc-pvqz_ecp.txt', 'r') as f:
                    for l in f.readlines()[2:]:
                        if '=' in l:
                            l_split = [s.strip() for s in l.split('=')]
                            entry.setdefault('E_all', {})[l_split[0]] = float(l_split[-1])
            except Exception:
                print('Error reading file ' + path_to_files + filename + 'wB97XD_ccsdt_cc-pvqz_ecp.txt')

        def read_Hcorr(entry, filename):
            try:
                with open(path_to_files + filename + 'wB97XD.td', 'r') as f:
                    entry.setdefault('E_all', {})['Hcorr'] = [float(l.split()[-1]) for l in f.readlines()
                                    if 'Thermal correction to Enthalpy' in l][0]
            except Exception:
                print('Error reading file ' + path_to_files + filename + 'wB97XD.td')

        for r in reference.reference:
            entry = {}
            entry['smi'] = reference.smiles[r]
            entry['name'] = reference.name[r]
            
            read_energies(entry, r)
            read_Hcorr(entry, r)

            entry['E_all']['dfH'] = reference.reference[r][0]
            entry['E_all']['u'] = reference.reference[r][1]

            self.ref[r] = entry

        entry = {}
        entry['smi'] = self.smi
        entry['name'] = self.name if self.name is not None else self.filename
        
        read_energies(entry, self.filename)
        read_Hcorr(entry, self.filename)

        self.ref[self.filename] = entry

    def init_gen(self):
        self.gen = ReactionGenerator(self.smi, [self.ref[name]['smi'] for name in self.ref 
                                                if name != self.filename], self.isodesmic)
        
    def run(self, max_rxns=1000, mp=False, n_proc=4):
        if mp:
            rxns = []
            
            with Pool(n_proc) as p:
                for res in p.imap_unordered(self.gen.generate, ((i, n_proc) for i in range(n_proc))):
                    rxns.extend(res)
        else:
            rxns = self.gen.generate(max_found=max_rxns)

        self.reactions = sorted(rxns, key=lambda x: sum(abs(x[1])))
        
        filenames = [self.filename] + [name for name in self.ref 
                                     if name != self.filename]

        with open(path_to_files + self.filename + 'wB97XD.all', 'w') as f:
            for rxn in self.reactions:
                entry = []
                entry.append(self.gen.reaction_to_str(*rxn))
                entry.append(sum(rxn[1]))
                entry.append(sum(abs(rxn[1])))
                entry.append(abs(entry[1]+entry[2]))

                E_sum = {}
                for r, c in zip(*rxn):
                    c_frac = c/rxn[1][0]
                    E_all = self.ref[filenames[r]]['E_all']

                    for n, e in E_all.items():
                        if n=='u':
                            E_sum[n] = E_sum.get(n, 0) + (c_frac*e)**2
                        else:
                            E_sum[n] = E_sum.get(n, 0) + c_frac*e
                
                E_core = E_sum['E(CV correction (TZ))']
                E_IT = E_sum['E(DLPNO-CCSD(T)/TZ-IT)']-E_sum['E(DLPNO-CCSD(T)/CC-PVTZ)']
                E_SO = E_sum['E(DKH correction)']

                E_r = E_sum['E(DLPNO-CCSD(T)/CBS)'] + E_SO + E_IT + E_core
                E_r_HF = E_sum['E(HF/CBS)'] + E_SO
                E_corr_r = E_r - E_r_HF
                
                H_r = E_r + E_sum['Hcorr']
                H_f = E_sum['dfH'] - H_r*ehtocal

                u = E_sum['u']**0.5

                entry.append(round(H_f, 3))
                entry.append(round(E_corr_r*ehtocal, 3))
                entry.append(round(E_sum['Hcorr']*ehtocal, 3))
                entry.append(entry[-1]+entry[-2])
                entry.append(round(u, 3))
                entry.append(round(2*u, 3))
                entry.append(round(H_r*ehtocal, 3))
                entry.append(0)
                entry.append(round(E_core*ehtocal, 3))
                entry.append(round(E_IT*ehtocal, 3))
                entry.append(0)
                entry.append(round(E_SO*ehtocal, 3))

                entry.append(self.gen.reaction_to_smi(*rxn))

                f.write(';'.join(map(str, entry)) + '\n')

if __name__ == '__main__':
    est = EnthalpyEstimator('m1_00_91-20-3_C10H8_', 'c1c2ccccc2ccc1', isodesmic=True)
    est.run(max_rxns=2000, mp=True)