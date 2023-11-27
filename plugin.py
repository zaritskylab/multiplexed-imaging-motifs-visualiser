from plugin_utils import *
import napari
import magicgui as mg
import imageio
import numpy as np
from matplotlib.backends.backend_agg import FigureCanvasAgg
import itertools


@mg.magicgui
class GraphVisualizer:
    def __init__(self):
        self.ui = self.create_ui()
        self.viewer = napari.Viewer()

    def create_ui(self):        
        self.csv_file = mg.widgets.FileEdit(filter="CSV files (*.csv)", mode="r", label="CSV File Path")
        self.patient_id_input = mg.widgets.LineEdit(name="Patient ID", label="Patient ID")
        self.fov_input = mg.widgets.LineEdit(name="FOV", label="Field of View (FOV)")
        self.visualize_button = mg.widgets.PushButton(text="Visualize")
        self.visualize_button.clicked.connect(self.on_button_click)
        self.image_file = mg.widgets.FileEdit(filter="Images (*.png *.xpm *.jpg *.bmp *.tif *.tiff)", mode="r", label="Image Path")
        self.motif_input = mg.widgets.Textarea(value="", label = "Motif Input")
        self.data_choice_widget = mg.widgets.ComboBox(choices=["onlyCSV", "segmentedImgae", "labeledImgae"], name="Data Type")
        layout = mg.widgets.Container(widgets=[self.csv_file, self.patient_id_input, self.fov_input, self.image_file, self.visualize_button, self.motif_input, self.data_choice_widget])
        return layout
    
    
    def on_button_click(self):
        print("Button clicked!")
        print(f"Motif input content: {self.motif_input.value}")
        csv_path = self.csv_file.value
        patient_id = self.patient_id_input.value
        fov = self.fov_input.value
        image_full_path = self.image_file.value
        if self.data_choice_widget.value == 'onlyCSV':
            data = pd.read_csv(fr'{csv_path}')
            if fov != None:
                data = data[(data['patient number'] == patient_id) & (data['fov'] == fov)].sample(frac = 0.05)
            else:
                data = data[(data['patient number'] == patient_id)]
            coords = data[['centroid-0', 'centroid-1']].values
            cell_types = data['pred'].values 
            tnbc_cells_type = generate_cell_type_structure(cell_types)
            print(len(np.unique(cell_types)))
            G_full, p2c,coords, coords_included = build_cell_graph(data, tnbc_cells_type, exclude=[])
           
        elif self.data_choice_widget.value == 'segmentedImgae':
            mapper = {'patints_col' : 'SampleID',
                        'cell_types_col' : 'cell_type',
                        'cell_index_col' : 'cellLabelInImage',
                        'fov_col' : 'FOV'}
            tb = pd.read_csv(fr'{csv_path}')
            cell_types = np.unique(tb['cell_type'])
            image_full = imageio.imread(image_full_path).astype(np.uint8)
            image_full = create_tagged_image(tb,image_full, mapper, patinet_number=patient_id, fov = fov)
            image_original, G_full, coords, point2cell_full, p2c = build_graph(image_full, list_of_cells_to_exclude = [4,5,6])
            tnbc_cells_type = generate_cell_type_structure(cell_types)
            print(tnbc_cells_type)
            
        elif self.data_choice_widget.value == 'labeledImgae':
            image_full = imageio.imread(image_full_path)  
            image_original, G_full, coords, point2cell_full, p2c = build_graph(image_full, list_of_cells_to_exclude = [4])
            tnbc = pd.read_csv(fr'{csv_path}')
            cell_types = np.unique(tnbc['CellType'])
            mapper_cells = {
            1: 'Unidentified',
            2: 'Endothelial',
            3: 'Mesenchyme',
            4: 'Tumor',
            5: 'Tregs',
            6: 'CD4 t cells',
            7: 'CD8 T cells',
            8: 'CD3 T cells',
            9: 'NK cells',
            10: 'B cells',
            11: 'Neutrophils',
            12: 'Macrophages',
            13: 'DC',
            14: 'DC/Mono',
            15: 'Mono/Neu',
            16: 'Immune other'
        }
            tnbc_cells_type = generate_cell_type_structure_from_tagged(mapper_cells)
            print(tnbc_cells_type)
        
        c = []
        ccc = []
        print(p2c)
        for i in p2c:
            cell_t = p2c[i]
            ccc.append(cell_t)
            c.append(tnbc_cells_type[cell_t]['color'])
        rotated_all_nodes_coords = rotate_point(coords, 270)
        self.viewer.add_points(rotated_all_nodes_coords, size=15, face_color=c, name="All Cells")
        node_ft_full = list(G_full.nodes(data=True))
        edges_full = list(G_full.edges())
        corrected_coords = []
        colors = []

        node_indices = [node[0] for node in node_ft_full]  # Extract the node indices
        # Only consider nodes that are part of the network subset
        for node, data in G_full.nodes(data=True):
            if node in node_indices:  # Check if node is part of the subset
                corrected_coords.append(coords[node])
                colors.append(tnbc_cells_type[data['cell_type']]['color'])

        
        rotated_coords = rotate_point(corrected_coords, 270)
        self.viewer.add_points(rotated_coords, size=15, face_color=colors, name="Only Commons")
        
        # Visualize edges as lines in Napari
        lines = []
        for edge in edges_full:
            lines.append([coords[edge[0]], coords[edge[1]]])
        
        rotated_lines = [rotate_point(line, 270) for line in lines]
        self.viewer.add_shapes(rotated_lines, shape_type='line', edge_color='white', name="Edges")
        #motifs = find_motif_in_graph(G_full, connections, non_connections)
        if not self.motif_input.value.strip():
            motif = {
                ('A', 'B'): {}, 
                ('A', 'C'): {},
                ('B', 'C'): {},
                'A': {'type': 9},
                'B': {'type': 11},
                'C': {'type': 11}
            }
        else:
            motif = parse_motif_input(self.motif_input.value)

        motifs = find_motifs(G_full, motif)
        all_motif_points = []
        all_motif_lines = []
        for motif_instance in motifs:
            # Extract the actual node indices from the motif instance
            motif_nodes = list(motif_instance.keys())
            motif_edges = list(itertools.combinations(motif_nodes, 2))
            # Rotate and visualize motif nodes
            rotated_motif_coords = [point for node in motif_nodes for point in rotate_point(coords[node], 270)]
            all_motif_points.extend(rotated_motif_coords)
            # Rotate and visualize motif edges
            rotated_motif_line_coords = []
            for line in motif_edges:
                if G_full.has_edge(line[0], line[1]):
                    start_point = rotate_point(coords[line[0]], 270)
                    end_point = rotate_point(coords[line[1]], 270)
                    rotated_motif_line_coords.extend([start_point, end_point])
            all_motif_lines.extend(rotated_motif_line_coords)
        print(len(all_motif_points))
        if len(all_motif_points) > 0:
            c = np.array(all_motif_points).reshape(-1,2)
            l = np.array(all_motif_lines).reshape(-1,2,2)
            self.viewer.add_shapes(l, shape_type='line', edge_color='red', name="Motif Edges")
            self.viewer.add_points(c, size=7, face_color='red', name="Motif Nodes", edge_color='red')

        ellipses = []
        colors = []
        labels = []
        y_coord = 30  # Starting y-coordinate for the center of the ellipse
        radius = 20  # Radius of each ellipse for both x and y axes
        for idx, details in tnbc_cells_type.items():
            ellipse_center = [y_coord, 30]
            ellipses.append([ellipse_center, [radius, radius]])
            
            colors.append(details['color'])
            labels.append(details['name'])
            y_coord += 60 

        # Add the Shapes layer with rectangles and colors
        shapes_layer = self.viewer.add_shapes(data=ellipses, shape_type='ellipse', edge_color=colors, face_color=colors, name="Legend")

        # Set the labels as the properties for the Shapes layer
        shapes_layer.properties = {'label': labels}
        shapes_layer.text = {'text': '{label}', 'size': 5, 'color': 'white', 'translation': [0, +7.5]}




visualizer = GraphVisualizer()


with napari.gui_qt():
    #viewer = napari.Viewer()
    visualizer.viewer.window.add_dock_widget(visualizer.ui, area='right')
