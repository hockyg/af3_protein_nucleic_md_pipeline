import os
import sys

def get_energies(sim1, sim3_start_idx, sim3_size, sim2, sim3):
    from openmm import unit
    c1 = sim1.context.getState(getEnergy=True,getPositions=True)
    total_energy = c1.getPotentialEnergy()
    system_positions = c1.getPositions()

    protein_positions = system_positions[:sim3_start_idx]
    nuc_positions = system_positions[sim3_start_idx:sim3_start_idx+sim3_size]

    sim2.context.setPositions(protein_positions)
    sim3.context.setPositions(nuc_positions)

    protein_energy = sim2.context.getState(getEnergy=True).getPotentialEnergy()
    nuc_energy = sim3.context.getState(getEnergy=True).getPotentialEnergy()
    #return energies in kJ/mol
    return [total_energy.value_in_unit(unit.kilojoule/unit.mole), protein_energy.value_in_unit(unit.kilojoule/unit.mole), nuc_energy.value_in_unit(unit.kilojoule/unit.mole), (total_energy-protein_energy-nuc_energy).value_in_unit(unit.kilojoule/unit.mole)]

def run_md(cif_file):
    def analyze(xtc, pdb):
        import mdtraj as md
        trj = md.load(xtc,top=pdb)
        #align to chain B, which is the big protein
        trj = trj.superpose(trj,atom_indices=trj.top.select('backbone and chainid 1'))
        #compute final RMSD and RMSF for chain A (designed molecule)
        rmsd = md.rmsd(trj,trj,atom_indices=trj.top.select('backbone and chainid 0'))[-1]
        rmsf = md.rmsf(trj,None,atom_indices=trj.top.select('backbone and chainid 0')).mean()
        return rmsd, rmsf

    def clean_pdb(input_pdb, output_pdb, keepMg=False):
        from openmm import app
        from pdbfixer import PDBFixer
        from openmm.app.element import magnesium
        fixer = PDBFixer(filename=input_pdb)
        chains = fixer.topology.chains()

        deleteAtoms = []
        MgPositions = []

        #this part is to remove 5' phosphorylation on nucleic acids
        for chainidx, chain in enumerate(chains):
            residues = [r for r in chain.residues()]
            #print(residues[0].name)
            if residues[0].name in ["A","C","T","G","U"]:
                atoms0 = residues[0].atoms()
                for atom in atoms0:
                    if atom.name in ["OP1","OP2","OP3","P"]:
                        deleteAtoms.append(atom)
            if residues[0].name in ["MG"]:
                atoms0 = residues[0].atoms()
                mg_idx = list(atoms0)[0].index
                MgPositions.append(fixer.positions[mg_idx])


        # delete these atoms

        modeller = app.Modeller(fixer.topology, fixer.positions)
        modeller.delete(deleteAtoms)

        fixer.topology = modeller.topology
        fixer.positions = modeller.positions
        fixer.removeHeterogens(False)

        if keepMg is True:
            # add back in magnesium in single chain with correct positions
            modeller = app.Modeller(fixer.topology, fixer.positions)
            c = modeller.topology.addChain("ION")
            for mgidx in range(len(MgPositions)):
                r = modeller.topology.addResidue("MG",c)
                modeller.topology.addAtom("MG",magnesium,r)
    
            for q in MgPositions:
                modeller.positions.append(q)

        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(7.0)

        app.PDBFile.writeFile(fixer.topology, fixer.positions, open(output_pdb, 'w'))
    
    import openmm
    from openmm import app
    from openmm import unit

    target=cif_file

    target_clean=target.replace('.cif','.autoclean.pdb')
    prefix=target.replace('.cif','_md')
    
    clean_pdb(target,target_clean)
    
    if not os.path.exists(f'{prefix}_run.xtc'):
        temperature=300
        #for production:
        equilibrationSteps=25000 #100 ps
        nsteps=2500000 # 10 ns

        #for testing
        #equilibrationSteps=2500 #10 ps
        #nsteps=25000 # 100ps 

        pair_energy_steps=1000
    
        forcefield = app.ForceField('amber14-all.xml', 'implicit/gbn2.xml')

        pdb = app.PDBFile(target_clean)

        modeller = app.Modeller(pdb.topology, pdb.positions)
        
        print("Setting up platform")
        platform = openmm.Platform.getPlatformByName('CUDA')
        #platform = openmm.Platform.getPlatformByName('CPU')
        #platform.precision='single'
        
        properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
        #properties = {}
        
        system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.NoCutoff,
                nonbondedCutoff=1*unit.nanometer, constraints=app.AllBonds, hydrogenMass=3*unit.amu)
        integrator = openmm.LangevinMiddleIntegrator(temperature*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
        simulation = app.Simulation(modeller.topology, system, integrator, platform, properties)
        simulation.context.setPositions(modeller.positions)


        #make other systems to calculat eenergy

        protein_modeller = app.Modeller(modeller.topology, modeller.positions)
        nuc_modeller = app.Modeller(modeller.topology, modeller.positions)

        not_protein_chains = [c for cidx, c in enumerate(protein_modeller.topology.chains()) if cidx!=0]
        not_nuc_chains = [c for cidx, c in enumerate(nuc_modeller.topology.chains()) if cidx!=1]

        protein_modeller.delete(not_protein_chains)
        nuc_modeller.delete(not_nuc_chains)

        #assume first protein then nuc
        n_protein_atoms = len([a for a in protein_modeller.topology.atoms()])
        n_nuc_atoms = len([a for a in nuc_modeller.topology.atoms()])

        systemProt = forcefield.createSystem(protein_modeller.topology, nonbondedMethod=app.NoCutoff,
                nonbondedCutoff=1*unit.nanometer, constraints=app.AllBonds, hydrogenMass=3*unit.amu)
        pintegrator = openmm.LangevinMiddleIntegrator(temperature*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
        simulationProt = app.Simulation(protein_modeller.topology, systemProt, pintegrator, platform, properties)
        simulationProt.context.setPositions(protein_modeller.positions)

        systemNuc = forcefield.createSystem(nuc_modeller.topology, nonbondedMethod=app.NoCutoff,
                nonbondedCutoff=1*unit.nanometer, constraints=app.AllBonds, hydrogenMass=3*unit.amu)
        nintegrator = openmm.LangevinMiddleIntegrator(temperature*unit.kelvin, 1/unit.picosecond, 0.004*unit.picoseconds)
        simulationNuc = app.Simulation(nuc_modeller.topology, systemNuc, nintegrator, platform, properties)
        simulationNuc.context.setPositions(nuc_modeller.positions)

        #simulation.reporters.append(app.PDBReporter(f'{prefix}_min.pdb', ,enforcePeriodicBox=False))

        print("Minimizing")
        simulation.minimizeEnergy()
        
        print('Equilibrating...')
        simulation.reporters.append(app.StateDataReporter(f'{prefix}.thermo.log', 1000, step=True,
                potentialEnergy=True, temperature=True, speed=True))
        simulation.reporters.append(app.PDBReporter(f'{prefix}_equil.pdb', equilibrationSteps ,enforcePeriodicBox=False))
        
        simulation.context.setVelocitiesToTemperature(temperature)
        simulation.step(equilibrationSteps)

        interaction_fh = open(f'{prefix}.interactions.log','w')
        interaction_fh.write("# Step P+N P N dE (kJ/mol)\n")
        
        #print("\tPost equilibration energies (P+N,P,N) (kJ/mol)")
        energy_list = []
        energies = get_energies(simulation, n_protein_atoms,  n_nuc_atoms, simulationProt, simulationNuc)
        interaction_fh.write("0 %f %f %f %f\n"%tuple(energies))
        interaction_fh.flush()

        print("Running MD")
        simulation.reporters.pop()
        xtcReporter = app.XTCReporter(f'{prefix}_run.xtc', pair_energy_steps, enforcePeriodicBox=False)
        simulation.reporters.append(xtcReporter)
        current_step = 0
        while current_step < nsteps:
            simulation.step(pair_energy_steps)
            energies = get_energies(simulation, n_protein_atoms,  n_nuc_atoms, simulationProt, simulationNuc)
            current_step = current_step + pair_energy_steps
            interaction_fh.write("%d %f %f %f %f\n"%tuple( [current_step]+energies ))
            interaction_fh.flush()
            energy_list.append(energies)

        interaction_fh.close()

#    rmsd, rmsf = analyze(f'{prefix}_run.xtc', f'{prefix}_equil.pdb')
    
#    print (rmsd, rmsf)
    return 0


def main():
    if len(sys.argv)<2:
        print(f"Usage: {sys.argv[0]} cif")
        sys.exit(1)

    cif_file = sys.argv[1]
    if not os.path.exists(cif_file):
        print(f"Error, {cif_file} does not exist")
        sys.exit(1)

    run_md(cif_file)

if __name__ == "__main__":
    main()
