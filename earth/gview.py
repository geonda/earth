"""
A.Geondzhian, N. Davydov
"""
import logging
import numpy as np

import plotly.graph_objects as go
import plotly.express as px

# from siman.core.structure import Structure


class visual:
    """
    class for crystal structure visualization
    """

    def __init__(self, atoms):
        """
        Initialize the visual
        Args:
        atoms  - for now only siman object which defines structure.
        """

        self.atoms = atoms
        self.param = dict(scatter=True,
                          projection=False,
                          projection_type='xy',
                          cell_vectors=False,
                          bonds=True,
                          bonds_color=False,
                          bonds_color_threshold=0.06,
                          bonds_elements=['Ni', 'O'],
                          bonds_length=True,
                          bonds_color_scale='Plasma',
                          bonds_threshold=2.4,
                          resolution=20,
                          )

        # self.fig = go.Figure()
        self.fig = None
        self.annotations = []
        self.bond_pairs = {}
        self.data = []
        self.abc = {0: "a", 1: "b", 2: 'c'}
        # test
        self.color_labels = {
            'Ni': {'color': 'rgb(77,0,242)', 'size': 0.3, 'opacity': 0.95},
            'O': {'color': 'rgb(235,0,0)', 'size': 0.15, 'opacity': 0.95},
            'Li': {'color': 'grey', 'size': 0.2, 'opacity': 0.5}
        }

    


    def plot(self, param=None):
        """Method to create and fill with content plotly figure. 

        Args:
            param (dict, optional): dictionary of optional parameters to adjust visualizatoin. 
            Defaults to None.
        """
        if param is not None:
            self.param.update(param)

        for index, (name, pos) in enumerate(zip(self.atoms.get_chemical_symbols(), self.atoms.positions)):
            self.data.append(self.create_atom_instance(index, name, pos))
        # print(self.data.append)data=self.data
        self.fig = go.Figure(data=self.data)

        self._add_bonds()
        self._add_cell()
        self._general_fomarting()
        self.fig.data = self.fig.data[::-1]

    def plot2d(self, index, name, pos):
        """method to draw atoms as circiles in 2d

        Args:
            index (int): index which will be reflected in the hoverinfo
            name (str): atomic symbol
            pos (list): 2d position of the atom

        Returns:
            plotly graphical object: instance for the given atom
        """
        return go.Scatter(x=[pos[0]],
                          y=[pos[1]],
                          opacity=self.color_labels[name]['opacity'],
                          name=f'{name} index {index}',
                          marker=dict(size=self.color_labels[name]['size']*100,
                                      color=self.color_labels[name]['color']),
                          showlegend=False)

    def plot3d(self, index, name, pos):
        """method to draw atoms spheres with _ms or points with go.Scatter3d

        Args:
            index (int): index which will be reflected in the hoverinfo
            name (str): atomic symbol
            pos (list): 2d position of the atom

        Returns:
            plotly graphical object: instance for the given atom
        """
        if self.param['scatter']:
            return go.Scatter3d(x=[pos[0]],
                                y=[pos[1]],
                                z=[pos[2]],
                                opacity=self.color_labels[name]['opacity'],
                                name=f'{name} index {index}',
                                marker=dict(size=self.color_labels[name]['size']*50,
                                            color=self.color_labels[name]['color']),
                                showlegend=False)
        else:
            (x_pns_surface, y_pns_surface, z_pns_suraface) = self._ms(
                pos, self.color_labels[name]['size'])
            return go.Surface(x=x_pns_surface, y=y_pns_surface, z=z_pns_suraface,
                              opacity=self.color_labels[name]['opacity'],
                              name=f'{name} index {index}',
                              colorscale=[[0, self.color_labels[name]['color']],
                                          [1, self.color_labels[name]['color']]],
                              showscale=False)

    def apply_projection(self, pos: list):
        """ Method to process projections.

        Args:
            index (int): index of a given atom
            name (str): atomic symbol
            pos (list): list of x,y,z coordinates in agnstrom 
        """
        if self.param['projection']:
            if 'xz' in self.param['projection_type'].lower():  # type: ignore
                projected_x, projected_y = pos[0], pos[2]
            elif 'yz' in self.param['projection_type'].lower():  # type: ignore
                projected_x, projected_y = pos[1], pos[2]
            elif 'xy' in self.param['projection_type'].lower():  # type: ignore
                projected_x, projected_y = pos[0], pos[1]
            else:
                logging.warning(
                    'Something went wrong in apply_projection, using "xy" instead.')
                projected_x, projected_y = pos[0], pos[1]
            return [projected_x, projected_y]
        else:
            return pos

    def create_atom_instance(self, index, name, pos):
        """wrapper on plot3d plot2d

        Args:
            index (int): index of a given atom
            name (str): atomic symbol
            pos (list): postion in angst for 2d or 3d cases
        """
        if not self.param['projection']:
            return self.plot3d(index, name, pos)
        else:
            pos2d = self.apply_projection(pos)
            return self.plot2d(index, name, pos2d)

    def get_bonds(self):
        """Method to get information about bond. Requires some adjustment.
        """
        for i, (name, position) in enumerate(zip(self.atoms.get_chemical_symbols(), self.atoms.positions)):
            if name in self.param['bonds_elements']:  # type: ignore
                for j, (second_name, second_position) in enumerate(zip(self.atoms.get_chemical_symbols()[i:], self.atoms.positions[i:])):
                    # type: ignore
                    if second_name in self.param['bonds_elements']:
                        bond_length = np.sqrt(
                            np.dot(position-second_position, position-second_position))
                        if bond_length < self.param['bonds_threshold'] and bond_length > 0.:
                            self.bond_pairs.update({f'{i}{j}':
                                                    {'atom1': name,
                                                     'atom2': second_name,
                                                     'bond_length': bond_length,
                                                     'begin':  position,
                                                     'end': second_position,
                                                     'color': 'green',
                                                     'annotation':
                                                     {'x': ((position+second_position)/2)[0],
                                                      'y': ((position+second_position)/2)[1],
                                                         'z': ((position+second_position)/2)[2],
                                                      }
                                                     }})

    def bond_color(self,):
        """ Method to get color scale and mark bonds with diffrenet lengths. 
        """
        scale = [self.bond_pairs[k]['bond_length'] for k in self.bond_pairs]
        self.bond_pairs = {k: v for k, v in sorted(
            self.bond_pairs.items(), key=lambda item: item[1]['bond_length'])}
        edge = min(scale)
        tmp = int((-min(scale)+max(scale))/self.param['bonds_color_threshold'])
        n_colors = tmp if tmp != 0 else 1
        if n_colors != 1:
            colors = px.colors.sample_colorscale(self.param['bonds_color_scale'], [
                                                 n/(n_colors - 1) for n in range(n_colors)])
            step = (-min(scale)+max(scale))/(n_colors-1)
        else:
            colors = ['green']
            step = max(scale)

        index = 0
        for k in self.bond_pairs.keys():
            # print(index)
            if self.bond_pairs[k]['bond_length'] < edge+step:
                self.bond_pairs[k]['color'] = colors[index]
            else:
                edge = edge+step
                index += 1
                self.bond_pairs[k]['color'] = colors[index]

    def _general_fomarting(self):
        """Options for general formatting of plotly figure.
        """
        self.fig.update_layout(coloraxis_showscale=False)

        if not self.param['projection']:
            self.fig.update_scenes(xaxis_visible=False,  # type: ignore
                                   yaxis_visible=False, zaxis_visible=False)
            self.fig.update_layout(scene=dict(bgcolor='white'),
                                   plot_bgcolor='rgb(0,0,0)',
                                   paper_bgcolor='rgb(0,0,0)',
                                   margin=dict(r=0, b=0, l=0, t=0),
                                   width=1000,)
            self.fig.update_layout(scene_camera=dict(up=dict(x=0, y=1., z=0),
                                                     eye=dict(x=0, y=0., z=0.5)),
                                   margin=dict(r=0, b=0, l=0, t=0),)
            self.fig.update_scenes(camera_projection_type='orthographic')
            self.fig.update_layout(scene=dict(
                annotations=self.annotations, aspectmode='data'))
        else:
            # self.fig.update_scenes(xaxis_visible=False,  # type: ignore
            #                        yaxis_visible=False, yaxis_showticklabels=False, xaxis_showticklabels=False)
            self.fig.update_layout(plot_bgcolor='rgb(0,0,0,0)',
                                   paper_bgcolor='rgb(0,0,0,0)',
                                   margin=dict(r=50, b=50, l=50, t=50),
                                   hovermode='closest',
                                   )
            self.fig.update_layout(yaxis_scaleanchor="x",)
            self.fig.update_xaxes(visible=False)
            self.fig.update_yaxes(visible=False)
            self.fig.update_layout(showlegend=False)
            # self.fig.update_layout(yaxis_visible=False, yaxis_showticklabels=False)
        self.fig.update_layout(coloraxis={'showscale': False})

    def _plot_bonds(self):
        """Internal method to draw bonds for 2d and 3d configurations.
        """
        for bond in self.bond_pairs:
            if not self.param['projection']:
                self.fig.add_trace(self._bond_line_3d(self.bond_pairs[bond]))
                if self.param['bonds_length']:
                    self.annotations.append(dict(
                        x=self.bond_pairs[bond]['annotation']['x'],
                        y=self.bond_pairs[bond]['annotation']['y'],
                        z=self.bond_pairs[bond]['annotation']['z'],
                        text=f"{np.round(self.bond_pairs[bond]['bond_length'],2)}",
                        showarrow=False))
            else:
                self.fig.add_trace(self._bond_line_2d(self.bond_pairs[bond]))
                if self.param['bonds_length']:
                    projected = self.apply_projection([self.bond_pairs[bond]['annotation']['x'],
                                                       self.bond_pairs[bond]['annotation']['y'],
                                                       self.bond_pairs[bond]['annotation']['z']])
                    self.fig.add_annotation(
                        x=projected[0],
                        y=projected[1],
                        text=f"{np.round(self.bond_pairs[bond]['bond_length'],2)}",
                        showarrow=False)

    def _add_bonds(self):
        """Internal method to add bonds on the figure. 
        """
        if self.param['bonds']:
            self.get_bonds()
            if self.param['bonds_color']:
                self.bond_color()
            self._plot_bonds()

    def _add_cell(self):
        """Internal method to add cell vecotors on the figure.
        """
        if self.param['cell_vectors']:
            for index, vector in enumerate(self.atoms.cell):
                if not self.param['projection']:
                    self.fig.add_trace(self._cell_vectors_3d(vector, index))
                else:
                    self.fig.add_trace(self._cell_vectors_2d(vector, index))

    def _cell_vectors_3d(self, vector, index):
        """Internal method to handle 3d unitcell vectors.

        Args:
            vector (list): end point for radius vector
            index (int): 0-a,1-b,2-c

        Returns:
            plotly.go: line and name of the vector
        """
        self.annotations.append(dict(x=vector[0]+0.5,
                                     y=vector[1]+0.5,
                                     z=vector[2]+0.5,
                                     text=f"{self.abc[index]}",
                                     showarrow=False))

        return go.Scatter3d(x=[0, vector[0]],
                            y=[0, vector[1]],
                            z=[0, vector[2]],
                            mode='lines',
                            name='cell',
                            line=dict(color='grey', width=5), showlegend=False
                            )

    def _cell_vectors_2d(self, vector, index):
        """Internal method to handle 3d unitcell vectors.

        Args:
            vector (list): end point for radius vector
            index (int): 0-a,1-b,2-c

        Returns:
            plotly.go: line and name of the vector
        """
        projected = self.apply_projection(vector)

        self.fig.add_annotation(x=projected[0]+0.2,
                                y=projected[1]+0.2,
                                text=f"{self.abc[index]}",
                                showarrow=False)

        return go.Scatter(
            x=[0, projected[0]],
            y=[0, projected[1]],
            mode='lines',
            name='cell',
            line=dict(color='grey', width=1), showlegend=False
        )

    def _bond_line_3d(self, bond_pair):
        """Internal method to draw bonds in 3d.

        Args:
            bond_pair (dict): instance of the bond_pair dict with info about this bond

        Returns:
            plotly.go : line instance
        """
        return go.Scatter3d(
            x=[bond_pair['begin'][0], bond_pair['end'][0]],
            y=[bond_pair['begin'][1], bond_pair['end'][1]],
            z=[bond_pair['begin'][2], bond_pair['end'][2]],
            mode='lines',
            opacity=0.4,
            name=f"{bond_pair['atom1']}-{bond_pair['atom2']}",
            text=f"{np.round(bond_pair['bond_length'],2)}",
            hoverinfo="text+name",
            line=dict(color=bond_pair['color'], width=7), showlegend=False
        )

    def _bond_line_2d(self, bond_pair):
        """Internal method to draw bonds in to 2d.

        Args:
            bond_pair (dict): instance of the bond_pair dict with info about this bond

        Returns:
            plotly.go : line instance
        """
        begin = self.apply_projection(bond_pair['begin'])
        end = self.apply_projection(bond_pair['end'])
        return go.Scatter(x=[begin[0], (begin[0]+end[0])/2, end[0]], y=[begin[1], (begin[1]+end[1])/2, end[1]],
                          mode='lines', opacity=0.4,
                          name=f"{bond_pair['atom1']}-{bond_pair['atom2']}",
                          text=f"{np.round(bond_pair['bond_length'],2)}",
                          hoverinfo="name+text",
                          line=dict(color=bond_pair['color'], width=7),
                          showlegend=False
                          )

    def _ms(self, pos, radius):
        """local method to draw atomic spheres

        Args:
            pos (list): postion of the center
            radius (float): radius of the sphere

        Returns:
            (x,y,z): mesh for the sphere
        """
        u, v = np.mgrid[0:2*np.pi:self.param['resolution']*2j,  # type: ignore
                        0:np.pi:self.param['resolution']*1j]  # type: ignore
        x = radius * np.cos(u)*np.sin(v) + pos[0]
        y = radius * np.sin(u)*np.sin(v) + pos[1]
        z = radius * np.cos(v) + pos[2]
        return (x, y, z)
