import os
import subprocess
import time
from earth.ocean2cif import ocean_util
import yaml

class input_processing:
    def __init__(self, filename=''):
        self.input_structure = filename
        self.ocean_preload=ocean_util(filename=self.input_structure)
        self.ocean_json=self.ocean_preload.main()
    
    def print_input(self):
        pretty_input=yaml.dump(self.ocean_json,default_flow_style=False)
        print(pretty_input)
    

