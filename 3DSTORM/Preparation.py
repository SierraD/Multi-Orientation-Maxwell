class file_preparation(object):
    """
    Sierra Dean
    RIKEN SPring-8 Center
    2024
    """
    
    def __init__(self, file_xy, file_xz, magnification=20, pixelsize_xy=230, pixelsize_xz=13, TS_dims = 2):
        """
        A technique to prepare the ThunderSTORM results obtained from a LSM Scan for 3D mapping and analysis.
        
        Attributes:
        file_xy & file_xz: str 
            The name of the ThunderSTORM results table obtained from the XY and XZ orientations. 
            Example) "XY_File.csv", "XZ_File.csv"
        magnification: int
            The magnification of the lens used with the LSM camera. If not specified, 20x magnification will be assumed. 
        pixelsize_xy: int 
            The pixel size of the XY data in nm, calculated using the camera's pixel size, divided by the magnification.
            Example) Hamamatsu Quest at 20x magnificaiton = 4600 [nm] / 20 = 230 
                    If not specified, this will be the value assumed.
        pixelsize_xz: int
            The Z step size in nm between images obtained from the LSM.
            If not specified, 13 nm will be assumed.
        TS_dims: 2 or 3
            The number of positional dimensions specified in ThunderSTORM using the Z-stage Offset Menu. 
            If the third dimension was previously specified with correct Z step, no voxel adjustments will be made.
            
        Return:
            None. Will modify the data established in place.
        """
        if (type(file_xy) or type(file_xz)) != str:
            raise TypeError("The file name of ThunderSTORM localizations must be specified as a string.")
        if (type(magnification)) != int:
            raise TypeError("The magnification must be specified as an integer.")
        if (type(pixelsize_xy) or type(pixelsize_xz)) != int:
            raise TypeError("The aspects of the voxel must be specified as an integer.")
        if TS_dims != 2 and TS_dims != 3:
            raise ValueError("The number of dimensions in the ThunderSTORM files should be 2 or 3.")
        self.name_xy = file_xy
        self.name_xz = file_xz
        self.xypix = pixelsize_xy
        self.zpix = pixelsize_xz
        self.magnification = magnification
        self.df_xy = pd.read_csv(self.name_xy)
        self.df_xz = pd.read_csv(self.name_xz)
        self.dfxy = pd.concat([self.df_xy["x [nm]"]*self.xypix,
                               self.df_xy["y [nm]"]*self.xypix,
                               self.df_xy["uncertainty [nm]"]*self.xypix,
                               self.df_xy["intensity [photon]"],
                               self.df_xy["sigma [nm]"]*self.xypix], 
                               keys=["X_XY", "Y_XY", "U_XY", "I_XY", "S_XY"], axis=1)
        self.dfxz = pd.concat([self.df_xz["x [nm]"]*self.xypix, 
                               self.df_xz["y [nm]"]*self.zpix,
                               self.df_xz["uncertainty [nm]"]*self.xypix, 
                               self.df_xz["uncertainty [nm]"]*self.zpix, 
                               self.df_xz["intensity [photon]"],
                               self.df_xz["sigma [nm]"]*self.xypix, 
                               self.df_xz["sigma [nm]"]*self.zpix], 
                               keys=["X_XZ", "Z_XZ", "U_X", "U_Z","I_XZ", "S_X", "S_Z"], axis=1)
        if TS_dims == 2:
            self.dfxy.insert(2, "Z_XY", self.df_xy["frame"]*self.zpix)
            self.dfxz.insert(1, "Y_XZ", self.df_xz["frame"]*self.xypix)
        elif TS_dims == 3:
            self.dfxy.insert(2, "Z_XY", self.df_xy["frame"])
            self.dfxz.insert(1, "Y_XZ", self.df_xz["frame"])
        self.xy_len = len(self.dfxy)
        self.xz_len = len(self.dfxz)
        return 
    
    def set_to_center(self):
        """
        A technique used to center the positional data around the origin.
        
        Return:
            None. Will modify the data established in place.
        """
        for i in self.dfxy.keys()[0:2]:
            self.dfxy[str(i)] = self.dfxy[str(i)]-(max(self.dfxy[str(i)]) + min(self.dfxy[str(i)]))/2
        for j in self.dfxz.keys()[0:2]:
            self.dfxz[str(j)] = self.dfxz[str(j)]-(max(self.dfxz[str(j)]) + min(self.dfxz[str(j)]))/2
        self.dfxy["Z_XY"] = self.dfxy["Z_XY"]-max(self.dfxy["Z_XY"])/2
        self.dfxz["Z_XZ"] = self.dfxz["Z_XZ"]-max(self.dfxz["Z_XZ"])/2
        return self
    
    def limiting(self, axis, orientation, limit, direction):
        """
        A technique used to limit the data to a specified positional region of interest.
        
        Attributes:
        axis: str "X", "Y", "Z"
            The axis along which the data will be limited.
        orientation: str "XY", "XZ"
            The data orientation which contains the data to be limited.
        limit: num
            The number value which is the limiting factor for the data.
        direction: str
            The direction in which the data will be limited. If "lesser" is specified, all values above the 
            limit will be removed, and vice versa if "greater" is specified.
            
        Return:
            None. Will modify the data established in place.
        """
        if (type(axis) or type(orientation)) != str:
            raise TypeError("The axis and orientation must be specified as a string.")
        if (type(limit) != float) and (type(limit) != int):
            raise TypeError("The limit must be specified as a number.")
        method = str(axis)+"_"+str(orientation)
        if orientation == "XY":
            if direction == "less" or "Less" or "lesser" or "Lesser":
                indexes = [idx for idx, value in enumerate(self.dfxy[method]) if value >= limit]
            elif direction == "more" or "More" or "greater" or "Greater":
                indexes = [idx for idx, value in enumerate(self.dfxy[method]) if value <= limit]
            self.dfxy = self.dfxy.drop(indexes).reset_index()
        if orientation == "XZ":
            if direction == "less" or "Less" or "lesser" or "Lesser":
                indexes = [idx for idx, value in enumerate(self.dfxz[method]) if value >= limit]
            elif direction == "more" or "More" or "greater" or "Greater":
                indexes = [idx for idx, value in enumerate(self.dfxz[method]) if value <= limit]
            self.dfxz = self.dfxz.drop(indexes).reset_index()
        return self
