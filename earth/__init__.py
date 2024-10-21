import os
import subprocess
import time
import matplotlib.pyplot as plt
import yaml
from earth.input_manager import input_processing
import ase.io
import spglib 
from ase.atoms import Atoms
from ase.visualize import view as avw
import earth.gview
import ase.build
import numpy as np
from earth.workflow_graph import workflow
from pathlib import Path
import json
from earth.ocean2cif import ocean_util


class default:
    def __init__(self):
        self.input={
        "dft": "qe",
        "dft_energy_range": 50,
        "ecut": "-1",
        "ecut.quality": "high",
        "opf.program": "hamann",
        "para_prefix": "mpirun -n 8",
        "pp_database": "ONCVPSP-PBE-PDv0.4-stringent",
        "screen_energy_range": 150,
        "cnbse.broaden": -1,
        "core_offset": "true",
        "cnbse.spect_range": "1000 -20 80",
        "screen.nkpt": "-2.55",
        "nkpt" : "-4",
        "ngkpt" : "-4",
        "screen.final.dr": "0.02",
        "screen.grid.rmax": "10",
        "screen.grid.rmode": "lagrange uniform",
        "screen.grid.ang": "5 11 11 9 7",
        "screen.grid.deltar": "0.10 0.15 0.25 0.25 0.25",
        "screen.grid.shells": " -1 4 6 8 10",
        "screen.lmax": "2",
        "cnbse.rad": "5.5",
        "screen.shells": "3.5 4.0 4.5 5.0 5.5 6.0",
        "cnbse.niter": 1000,
        "haydock_convergence": " 0.001 5 "
    }
    
class manage_structure:
    def __init__(self,filename=None,workflow=None):
        self.filename=filename
        self.makePrimitive=True
        self.atoms=ase.io.read(filename)
        self.ocean_atoms=self.convert_to_ocean()
        self.workflow=workflow
        self.workflow.add_node(str(self.atoms.symbols))
        self.Elements=self.ocean_atoms.symbols


    def convert_to_ocean(self):
        new_unit_cell, new_scaled_positions, new_numbers = spglib.standardize_cell((self.atoms.cell, self.atoms.get_scaled_positions(), self.atoms.numbers ), to_primitive=self.makePrimitive, symprec=5e-3 )
        ocean_atoms = Atoms(new_numbers, cell=new_unit_cell, scaled_positions=new_scaled_positions)
        return ocean_atoms
        
    def view(self, mode='gview'):
        if mode=='ase':
            avw(self.ocean_atoms)
        elif mode=='supercell-gview':
            self.supercell=sp=ase.build.make_supercell(self.atoms, np.eye(3)*3, wrap=True, order='cell-major', tol=1e-05)
            self.svs=gview.visual(self.supercell)
            self.svs.plot(param=dict( cell_vectors=True,
                                                projection=True,))
            self.svs.fig.show()
        else:
            self.vs=gview.visual(self.atoms)
            self.vs.plot(param=dict( cell_vectors=True,
                                        projection=True,))
            self.vs.fig.show()

    

    
        
class ocean_wrapper:
    def __init__(self,filename):
        self.path={}
        self.workflow=workflow()
        self.structure=manage_structure(filename=filename, workflow=self.workflow)
        self.path[list(self.workflow.graph.nodes)[0]]=os.getcwd()
        for item in set(self.structure.atoms.get_chemical_symbols()):
            self.path[item]=dict()
        
    def xas(self,element=None,edge=None):
        first_node=list(self.workflow.graph.nodes)[0]
        self.workflow.add_instance(first_node,element)
        self.workflow.add_instance(element,f'{element}-{edge} edge')
        self.workflow.graph.nodes(data=True)
        self.path[element].update({edge:os.path.join(os.getcwd(),f'{first_node}/{element}/{edge}/')})
        Path(f'{first_node}/{element}/{edge}').mkdir( parents=True, exist_ok=True)
        ocean_json=default().input
        ocean_json["edges"]=f"-{element} {edge}"
        ocean_util().write_ocean_in(filename=f"{self.path[element][edge]}/ocean.in",
                                   atoms=self.structure.atoms,
                                     input_data=ocean_json)
   