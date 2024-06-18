class preparation(object):
    """
    This file is part of the 3D STORM software
    
    File author(s): Sierra Dean <ccnd@live.com>
    
    Distributed under the GPLv3 Licence.
    See accompanying file LICENSE.txt or copy at
        http://www.gnu.org/licenses/gpl-3.0.html
        
    source: https://github.com/SierraD/3DSTORM
    
    Last Updated: June 18 2024
    """
    def __init__(self):
        """
        A technique to prepare multi two-dimensional ThunderSTORM analysis files for three-dimensional
        analysis. 
        
        This method requires two different orientation files obtained from ThunderSTORM, containing
        two-dimensional pixel-value localizations, and requires the user to input the voxel size
        of the images by defining the in-plane pixel size, as well as the in-depth step size between 
        subsequent images. 
        
        The values in the ThunderSTORM files will then be corrected from pixel size to the nanometer
        scale, and can be set to a zero position, or limited to a user-defined region of interest, 
        as well as downloaded as a CSV file.
        
        Attributes:
            None.
            
        Return: 
            None.
        """
        
        return 
    
    def setting(self, file_xy, file_xz, magnification=20, pixelsize_xy=230, pixelsize_xz=13, TS_dims = 2):
        """
        A technique to download two different orientation two-dimensional ThunderSTORM analysis files
        and correct the scale from pixel size to nanometers.
                
        Attributes:
        file_xy & file_xz: str "XY_File.csv", "XZ_File.csv", etc.
            The name of the ThunderSTORM results table obtained from the XY and XZ orientations. 
        magnification: int
            The magnification of the lens used to obtain the images. If not specified, 20x magnification 
            will be assumed. 
        pixelsize_xy: int 
            The pixel size of the XY data in nm, calculated using the camera's pixel size, divided by 
            the magnification.
            Example) Hamamatsu Quest at 20x magnificaiton = 4600 [nm] / 20 = 230 [nm] 
                    If not specified, this will be the value assumed.
        pixelsize_xz: int
            The Z step size in nm between images obtained from the LSM.
            If not specified, 13 nm will be assumed.
        TS_dims: 2 or 3
            The number of positional dimensions specified in ThunderSTORM using the Z-stage Offset Menu. 
            If the third dimension was previously specified with correct Z step, no voxel adjustments will 
            be made.
            
        Return:
            None. Will modify the data established in place.
        """
        self.name_xy = file_xy
        self.name_xz = file_xz
        self.xypix = pixelsize_xy
        self.zpix = pixelsize_xz
        self.magnification = magnification
        self.df_xy = pandas.read_csv(self.name_xy)
        self.df_xz = pandas.read_csv(self.name_xz)
        self.dfxy = pandas.concat([self.df_xy["x [nm]"]*self.xypix,
                               self.df_xy["y [nm]"]*self.xypix,
                               self.df_xy["uncertainty [nm]"]*self.xypix,
                               self.df_xy["intensity [photon]"],
                               self.df_xy["offset [photon]"],
                               self.df_xy["bkgstd [photon]"],
                               self.df_xy["sigma [nm]"]*self.xypix], 
                               keys=["X_XY", "Y_XY", "U_XY", "I_XY", "O_XY", "B_XY", "S_XY"], axis=1)
        self.dfxz = pandas.concat([self.df_xz["x [nm]"]*self.xypix, 
                               self.df_xz["y [nm]"]*self.zpix,
                               self.df_xz["uncertainty [nm]"]*self.xypix, 
                               self.df_xz["uncertainty [nm]"]*self.zpix, 
                               self.df_xz["intensity [photon]"],
                               self.df_xz["offset [photon]"],
                               self.df_xz["bkgstd [photon]"],
                               self.df_xz["sigma [nm]"]*self.xypix, 
                               self.df_xz["sigma [nm]"]*self.zpix], 
                               keys=["X_XZ", "Z_XZ", "U_X", "U_Z","I_XZ", "O_XZ", "B_XZ", "S_X", "S_Z"], axis=1)
        if TS_dims == 2:
            self.dfxy.insert(2, "Z_XY", (self.df_xy["frame"]*self.zpix))
            self.dfxz.insert(1, "Y_XZ", (self.df_xz["frame"]*self.xypix))
        elif TS_dims == 3:
            self.dfxy.insert(2, "Z_XY", self.df_xy["frame"])
            self.dfxz.insert(1, "Y_XZ", self.df_xz["frame"])
        self.xy_len = len(self.dfxy)
        self.xz_len = len(self.dfxz)
        return self
    
    def set_to_center(self):
        """
        A technique to center the data to a zero center position.
        
        Attributes:
            None.
        
        Return:
            None. Will modify the data established in place.
        """
        self.dfxz["X_XZ"] = self.dfxz["X_XZ"]-(max(self.dfxy["X_XY"])+min(self.dfxy["X_XY"]))/2
        self.dfxy["X_XY"] = self.dfxy["X_XY"]-(max(self.dfxy["X_XY"])+min(self.dfxy["X_XY"]))/2
        self.dfxz["Y_XZ"] = self.dfxz["Y_XZ"]-(max(self.dfxy["Y_XY"])+min(self.dfxy["Y_XY"]))/2
        self.dfxy["Y_XY"] = self.dfxy["Y_XY"]-(max(self.dfxy["Y_XY"])+min(self.dfxy["Y_XY"]))/2
        self.dfxy["Z_XY"] = self.dfxy["Z_XY"]-max(self.dfxz["Z_XZ"])/2
        self.dfxz["Z_XZ"] = self.dfxz["Z_XZ"]-max(self.dfxz["Z_XZ"])/2
        return self
    
    
    def limiting(self, axis, limit, direction):
        """
        A technique to limit the data to a specified positional region of interest.
        
        Attributes:
        axis: str "X", "Y", "Z"
            The axis along which the data will be limited.
        limit: num
            The number value which is the limiting factor for the data.
        direction: str "Less", "Lesser", "More", "Greater", etc.
            The direction in which the data will be limited. If "lesser" is specified, all values above the 
            limit will be removed, and vice versa if "greater" is specified.
            
        Return:
            None. Will modify the data established in place.
        """
        if direction == ("less" or "Less" or "lesser" or "Lesser"):
            indexes_xy = [idy for idy, valuey in enumerate(self.dfxy[axis+"_XY"]) if valuey >= limit]
            indexes_xz = [idz for idz, valuez in enumerate(self.dfxz[axis+"_XZ"]) if valuez >= limit]
            self.dfxy = self.dfxy.drop(indexes_xy).reset_index(drop=True)
            self.dfxz = self.dfxz.drop(indexes_xz).reset_index(drop=True)
        elif direction == ("more" or "More" or "greater" or "Greater"):
            indexes_xy = [idy for idy, valuey in enumerate(self.dfxy[axis+"_XY"]) if valuey <= limit]
            indexes_xz = [idz for idz, valuez in enumerate(self.dfxz[axis+"_XZ"]) if valuez <= limit]
            self.dfxy = self.dfxy.drop(indexes_xy).reset_index(drop=True)
            self.dfxz = self.dfxz.drop(indexes_xz).reset_index(drop=True)
        return self
    
    def download_dataframe(self, filename="3DSTORM"):
        """
        A technique to download the data prepared by the preparation.py method as a CSV file.
        
        Attributes:
            None.
        Return:
            None. Will download the dataframe as a CSV file.
        """
        download_df = pandas.concat([self.dfxy,self.dfxz],axis=1,sort=False)
        download_df = download_df.rename(columns={'X_XY': 'x_xy [nm]', 
                                                  'Y_XY': 'y_xy [nm]', 
                                                  'Z_XY': 'z_xy [nm]', 
                                                  'U_XY': 'uncertainty_xy [nm]', 
                                                  'I_XY': 'intensity_xy [photon]',
                                                  'O_XY': 'offset_xy [photon]',
                                                  'B_XY': 'bkgstd_xy [photon]',
                                                  'S_XY': 'sigma_xy [nm]', 
                                                  'X_XZ': 'x_xz [nm]', 
                                                  'Y_XZ': 'y_xz [nm]', 
                                                  'Z_XZ': 'z_xz [nm]', 
                                                  'U_X': 'uncertainty_x [nm]',
                                                  'U_Z': 'uncertainty_z [nm]',
                                                  'I_XZ': 'intensity_xz [photon]',
                                                  'O_XZ': 'offset_xz [photon]',
                                                  'B_XZ': 'bkgstd_xz [photon]',
                                                  'S_X': 'sigma_x [nm]', 
                                                  'S_Z': 'sigma_z [nm]'})
        download_df.to_csv(filename+".csv", index=False, encoding='utf-8')
        return self