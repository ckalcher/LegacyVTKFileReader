#### LegacyVTKFileReader ####
# For structured grid with point data

import enum
import re
from typing import Any

import numpy as np
from ovito.data import *
from ovito.io import FileReaderInterface
from ovito.vis import *
from traits.api import *


class LegacyVTKFileReader(FileReaderInterface):

    class GridType(enum.StrEnum):
        CELLDATA = "Cell data"
        POINTDATA = "Point data"

    grid_type = Enum(GridType.CELLDATA, [GridType.CELLDATA, GridType.POINTDATA], label="Voxel grid type")  

    @staticmethod
    def detect(filename: str):
        try:
            with open(filename, "r") as f:
                line = f.readline()
                return line.strip() == "# vtk DataFile Version 2.0"
        except OSError:
            return False

    def parse(self, data: DataCollection, filename: str, frame_info: Any, **kwargs: Any):
        with open(filename, "r") as f:

            # Initialize data storage
            dimensions = None
            origin = None
            spacing = None
            scalars = []
            num_points = None

            # Regular expressions for parsing
            dim_pattern = re.compile(r"DIMENSIONS (\d+) (\d+) (\d+)")
            origin_pattern = re.compile(r"ORIGIN ([\d\.\-eE]+) ([\d\.\-eE]+) ([\d\.\-eE]+)")
            spacing_pattern = re.compile(r"SPACING ([\d\.\-eE]+) ([\d\.\-eE]+) ([\d\.\-eE]+)")
            scalars_pattern = re.compile(r"SCALARS (\w+) (\w+)")
            num_points_pattern = re.compile(r"POINT DATA (\d+)")

            # Read and parse the file
            with open(filename, 'r') as file:
                lines = file.readlines()
                for line in lines:
                    # Extract dimensions
                    dim_match = dim_pattern.match(line)
                    if dim_match:
                        dimensions = tuple(map(int, dim_match.groups()))
                        nx, ny, nz = dimensions
                    # Extract origin
                    origin_match = origin_pattern.match(line)
                    if origin_match:
                        origin = tuple(map(float, origin_match.groups()))

                    # Extract spacing
                    spacing_match = spacing_pattern.match(line)
                    if spacing_match:
                        spacing = tuple(map(float, spacing_match.groups()))

                    num_points_match = num_points_pattern.match(line)
                    if num_points_match:
                        num_points = map(int, num_points_match.groups())

                    # Start reading scalar values after encountering 'LOOKUP_TABLE'
                    if line.strip() == "LOOKUP_TABLE default":
                        # Scalar values follow after the LOOKUP_TABLE line
                        scalar_lines = lines[lines.index(line) + 1:]
                        for scalar_line in scalar_lines:
                            values = scalar_line.split()
                            if values:
                                scalars.extend(map(float, values))

            # Convert scalar values into a numpy array
            field_data = np.array(scalars)
            #assert len(field_data) == num_points

            # Optionally reshape the scalars array to match the grid dimensions
            if dimensions:
                field_data = field_data.reshape(dimensions)
           
            # Create a new SimulationCell object defining the outer spatial dimensions
            # of the grid and the boundary conditions, and add it to the DataCollection:
            cell = data.create_cell(
                matrix=[[nx*spacing[0],0,0,origin[0]],[0,ny*spacing[1],0,origin[1]],[0,0,nz*spacing[2],origin[2]]],
                pbc=(True, True, True)
            )

            # Create the VoxelGrid object
            grid = data.grids.create(
                identifier="density_field",
                domain=cell,
                shape=(nx,ny,nz),
                vis=VoxelGridVis(enabled=True, transparency=0.6)
            )
            if self.grid_type == LegacyVTKFileReader.GridType.POINTDATA:
                grid.grid_type = VoxelGrid.GridType.PointData
                
            grid.create_property('Field Value', data=field_data.flatten(order='F'))
