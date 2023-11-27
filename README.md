# Motif Visualizer for Multiplexed Imaging

## Introduction

**Motif Visualizer for Multiplexed Imaging** is a specialized tool designed to enhance the analysis and visualization of multiplexed imaging data. Multiplexed imaging is a cutting-edge technique in the field of biomedical research and diagnostics, enabling the simultaneous detection and visualization of multiple biomarkers within a single sample. This technique is pivotal in understanding complex biological systems and diseases by allowing researchers to observe the spatial distribution and interaction of multiple biomolecules in tissues or cells.

This tool specifically focuses on spatial network motifs within multiplexed imaging data. Spatial network motifs are recurring, significant patterns of spatial relationships between different types of biomolecules. Identifying these motifs is crucial for understanding the cellular architecture and its functional implications in health and disease.

## Features

- **Interactive View:** A dynamic, user-interactive visualization interface that allows for real-time exploration and analysis of spatial network motifs in multiplexed imaging data.
- **Static View:** A high-resolution, static image output that captures the intricacies of spatial relationships in a comprehensive and detailed manner.

## Installation and Usage

To get started with the Motif Visualizer for Multiplexed Imaging, follow these steps:

1. Clone the repository:
   ```
   git clone [repository URL]
   ```
2. Navigate to the cloned directory and install the required dependencies:
   ```
   cd [repository directory]
   pip install -r requirements.txt
   ```
3. Run the visualizer:
   ```
   python visualizer.py
   ```

Replace `[repository URL]` and `[repository directory]` with your actual repository URL and directory name.

## Examples

Here are some examples of how the Motif Visualizer can be used to explore and understand multiplexed imaging data:

1. **Interactive View Example:**
   <p float="left">
  <img src="https://github.com/YuvalTamir2/multiplexed-imaging-motifs-visualiser/assets/72014577/aec5f115-e533-4d5b-ae6d-4e1d6e12ae58" width="900" height = "450" /> 
  <img src="https://github.com/YuvalTamir2/multiplexed-imaging-motifs-visualiser/assets/72014577/7531bfa2-32fc-4383-865d-038282a5e944" width="350" height = "350" />
   </p>
   ![Interactive View](path/to/interactive_view_image.png)

2. **Static View Example:**
   ```
   python static_view.py path/to/cells.csv path/to/segmented_image.tiff patient_number 
   ```
   ![Static View]![graph_overlay_static](https://github.com/YuvalTamir2/multiplexed-imaging-motifs-visualiser/assets/72014577/e649f0b5-53da-4a80-bbf9-cb0b8edd15ee)


Replace `path/to/interactive_view_image.png` and `path/to/static_view_image.png` with the paths to your example images.

## Contributing

Contributions to the Motif Visualizer for Multiplexed Imaging are welcome. Please read our contribution guidelines for more information.

## License

This project is licensed under [License Name] - see the LICENSE file for details.
