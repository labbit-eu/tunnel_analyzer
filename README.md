# Tunnel Analyzer

## Overview
The Tunnel Analyzer project is designed to analyze 3D Voronoi diagrams and quasi triangulations in QTF file format of protein snapshots obtained from molecular dynamics simulations. It consists of a C++ program responsible for refining the Voronoi diagram and conducting a search for paths between user-specified or automatically detected starting points on the protein surface. Additionally, it includes a Python script for path clustering, which categorizes paths from all snapshots into similar groups.

## Features
- **C++ Program**:
  - Analyzes and refines 3D Voronoi diagrams.
  - Conducts exhaustive search for paths between specified or detected starting points on protein surface.
  - Supports multiple input QTF files.
  
- **Python Path Clustering Script**:
  - Categorizes paths from all snapshots into similar groups.
  - Enables visual analysis of path variations and similarities across different protein snapshots.

## Requirements
- **C++ Compiler**: Required to compile and run the C++ program.
  - **CGAL Library** 
- **Python**: Required to execute the path clustering script.
  - **NumPy**
  - **scikit-learn**
- **QTF Files**: Input files containing 3D Voronoi diagrams and quasi triangulations of protein snapshots.
- **input.txt**: Input file containing configuration. Exemplary input file is in the repository. First line is a path to a directory containing QTF files. Second line is minimal radius used in Voronoi diagram refinement. Following lines are optional and correspond to user-defined starting points.  

## Output
- The C++ program outputs sets of paths from each QTF file in single JSON file - output.json
- The Python script categorizes these paths into similar groups. Each cluster is saved to a seperate PDB file.
