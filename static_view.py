import updated_napari_ui
import pandas as pd
import numpy as np
import imageio
import matplotlib.cm as cm

class tissueData():
    def __init__(self, cells_data_csv, cells_data_image_path, mapper):
        self.cell_data = cells_data_csv
        self.image = imageio.imread(cells_data_image_path)
        
        self.mapper_dict = mapper


    def create_tagged_image(self, patinet_number, fov = None):
        # Load the CSV data
        if fov is None:
            cell_data_specific_p = self.cell_data[self.cell_data[self.mapper_dict['patints_col']] == patinet_number]
        else:
            cell_data_specific_p = self.cell_data[(self.cell_data[self.mapper_dict['patints_col']] == patinet_number) \
                                                  & (self.cell_data[self.mapper_dict['fov_col']] == fov)]

        segmented_image_array = np.array(self.image)
        # Create a mapping from cell type to integer
        cell_types = cell_data_specific_p[self.mapper_dict['cell_types_col']].unique()
        cell_type_to_int = {cell_type: idx + 1 for idx, cell_type in enumerate(cell_types)}

        # Create a new image array for cell types
        cell_type_image_array = np.zeros_like(segmented_image_array)
        # Populate the new image array based on cell type mapping
        for _, row in cell_data_specific_p.iterrows():
            cell_label = row[self.mapper_dict['cell_index_col']]
            cell_type_int = cell_type_to_int[row[self.mapper_dict['cell_types_col']]]
            cell_type_image_array[segmented_image_array == cell_label] = cell_type_int
        return cell_type_image_array.astype(np.uint8)
    
    def generate_cell_type_structure(self):
        """
        Generate a dictionary structure similar to tnbc_cells_type for the provided cell types.
        """
        cell_types = self.cell_data[self.mapper_dict['cell_types_col']]
        unique_cell_types = np.unique(cell_types)
        num_unique_types = len(unique_cell_types)
        
        # Generate a wide range of unique colors using a colormap
        colormap = cm.get_cmap('tab20c', num_unique_types)
        colors = [colormap(i) for i in range(num_unique_types)]
        
        # Background and Unknown are hardcoded
        cell_type_structure = {
            #0: {'name': 'Background', 'color': 'black'},
            #1: {'name': 'Unknown', 'color': 'black'}
        }
        
        # Populate the dictionary with unique cell types and their colors
        for i, cell_type in enumerate(unique_cell_types):#, start=2):
            cell_type_structure[i] = {'name': cell_type, 'color': colors[i]}
            
        return cell_type_structure
    
    def generate_cell_type_structure_from_tagged(self, mapper):
        num_unique_types = len(mapper)
        
        colormap = cm.get_cmap('tab20c', num_unique_types)
        colors = [colormap(i) for i in range(num_unique_types)]
        
        # Background and Unknown are hardcoded
        cell_type_structure = {
            #0: {'name': 'Background', 'color': 'black'},
            #1: {'name': 'Unknown', 'color': 'black'}
        }
        
        # Populate the dictionary with unique cell types and their colors
        for i, cell_type in enumerate(mapper.keys()):#, start=2):
            cell_type_structure[cell_type] = {'name': mapper[cell_type], 'color': colors[i]}
            
        return cell_type_structure
    
    def main(self,
            patinet_number,
            fov,
            motif,
            is_tagged = False,
            list_of_cells_to_exclude = []):
        

        if is_tagged == False:
            self.image = self.create_tagged_image(patinet_number, fov)
            cell_types_dict = self.generate_cell_type_structure()
        else:
            cell_types_dict = self.generate_cell_type_structure_from_tagged(is_tagged)
        print(cell_types_dict)
        self.image_original, self.G_full, self.coords_full, self.point2cell_full, p2c = updated_napari_ui.build_graph(self.image, list_of_cells_to_exclude)
        motifs = updated_napari_ui.find_motifs(self.G_full, motif = motif)
        fig = updated_napari_ui.vis_graph_and_motifs(coords_full = self.coords_full,
                                                    subgraphs_8_list =  motifs,
                                                    image_original = self.image_original,
                                                    G_full = self.G_full,
                                                    point2cell_full = self.point2cell_full,
                                                    cell_types=cell_types_dict)
        


#########################################################################################################
#### params specific for the dataset 

# mapping col names
mapper = {'patints_col' : 'SampleID',
          'cell_types_col' : 'CellType',
          'cell_index_col' : 'cellLabelInImage',
          'fov_col' : 'FOV'}
# the motif we want to find
motif = {
    ('A', 'B'): {}, 
    ('A', 'C'): {},
    ('B', 'C'): {},
    'A': {'type': 9},
    'B': {'type': 10},
    'C': {'type': 10}
}
# if the image is tagged, provide the index:cell_type mapper
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

cells_data = pd.read_csv((r'../TNBC_morph/ONLY_CELLS.csv'))

visualizer = tissueData(cells_data, 
                        r"CellTypes.tiff",
                        mapper)

visualizer.main(patinet_number=1, 
                fov=None,
                motif=motif,
                is_tagged=mapper_cells,
                list_of_cells_to_exclude=[i for i in range(1,17)]
                )


