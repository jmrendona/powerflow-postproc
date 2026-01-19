import vtk
from vtk.util import numpy_support
import numpy as np
import os

class VTKExporter:
    
    def __init__(self, Xg, Yg, field, fieldname: str):
        
        self.fieldname = fieldname
        self.Xg = Xg
        self.Yg = Yg
        self.field = field
        
        self.nx, self.ny = Xg.shape
        self.dx = Xg[0,1] - Xg[0,0]
        self.dy = Yg[1,0] - Yg[0,0]
        
    def SaveVTI(self,filename: str):
        
        image = vtk.vtkImageData()
        image.SetDimensions(self.nx, self.ny, 1)
        image.SetSpacing(self.dx, self.dy, 1)
        image.SetOrigin(0,0,0)
        
        vtk_data = numpy_support.numpy_to_vtk(self.field.ravel(order='C'), deep=True)
        vtk_data.SetName(self.fieldname)
        
        image.GetPointData().SetScalars(vtk_data)
        
        writer = vtk.vtkXMLImageDataWriter()
        writer.SetFileName(filename + '.vti')
        writer.SetInputData(image)
        writer.Write()
        
        print(f'VTK file saved as {filename}.vti')
