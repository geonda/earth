import plotly.graph_objects as go
import networkx as nx

class workflow:
    def __init__(self):
        self.graph = nx.Graph()
        self.id=0

    def add_node(self, node):
        self.graph.add_node(node,data='test')

    def add_edge(self, node1, node2):
        self.graph.add_edge(node1, node2)

    def add_instance(self,node1,node2):
        self.id+=1
        self.graph.add_node(node2,data=f'test{self.id}')
        self.add_edge(node1,node2)
        

    def show(self):
        """Visualize the graph using Plotly."""
        pos = nx.shell_layout(self.graph)
        
        # Create edge traces
        edge_x = []
        edge_y = []
        
        for edge in self.graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            edge_x.append(x0)
            edge_x.append(x1)
            edge_x.append(None)  # None separates edges
            edge_y.append(y0)
            edge_y.append(y1)
            edge_y.append(None)  # None separates edges
        
        # Create node traces
        node_x = []
        node_y = []
        
        for node in self.graph.nodes():
            x, y = pos[node]
            node_x.append(x)
            node_y.append(y)
        annotations = []
    
        for edge in self.graph.edges():
            x0, y0 = pos[edge[0]]
            x1, y1 = pos[edge[1]]
            
            annotations.append(dict(
                ax=x0, ay=y0,
                axref='x', ayref='y',
                x=x1, y=y1,
                xref='x', yref='y',
                showarrow=True,
                arrowhead=2,  # Arrowhead style
                arrowsize=5,
                arrowcolor='gray'
            ))

        # Create figure
        fig = go.Figure()
        
        # Add edges to figure
        fig.add_trace(go.Scatter(x=edge_x, y=edge_y,
                                 line=dict(width=0.5, color='#888'),
                                 mode='lines'))
        
        
        # Add nodes to figure
        fig.add_trace(go.Scatter(x=node_x, y=node_y,
                                 mode='markers+text',
                                 text=list(self.graph.nodes()),
                                 textposition="top center",
                                 hoverinfo='text',
                                 marker=dict(
                                             size=30,
                                             color='#6175c1',
                                             line_width=2)))

        fig.update_layout(showlegend=False,
                          hovermode='closest',
                          annotations=annotations,
                          margin=dict(b=0,l=0,r=0,t=0),
                          xaxis=dict(showgrid=False,
                                     zeroline=False,
                                     showticklabels=False),
                          yaxis=dict(showgrid=False,
                                     zeroline=False,
                                     showticklabels=False))

        fig.show()


