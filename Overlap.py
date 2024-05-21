class overlap(object):
    """
    This file is part of the 3D STORM software
    
    File author(s): Sierra Dean <ccnd@live.com>
    
    Distributed under the GPLv3 Licence.
    See accompanying file LICENSE.txt or copy at
        http://www.gnu.org/licenses/gpl-3.0.html
        
    source: https://github.com/SierraD/3DSTORM
    
    Last Updated: May 21 2024
    """
    
    def __init__(self, data):
        """
        A technique to use two different orientations of ThunderSTORM results to determine the overlapping 
        points in 3D space.
        
        This method splits the calculation into determining the indexes of the overlap, then determining
        the localizations using the indexes, which is computationally easier than performing both 
        within the same method.
        
        Attributes:
        data: 
            The data previously developed and contained within the Preparations.py class.
            
        Return:
            None. Will modify the data established in place.
        """
        self.data = data
        self.zpix = self.data.zpix
        self.xypix = self.data.xypix
        self.magnification = self.data.magnification
        return 
    
    def indexes(self, z_range, xy_range):
        """
        A technique to determine the indexes of the dataframe where the points from the two different 
        orientations overlap in 3D space. 
        
        For each localization, a 3D sphere is defined around each position, with the radius in each direction
        determined by that point's uncertainty value in that direction. For uncertain values (i.e. Z in XY and 
        Y in XZ), the uncertainty is taken as the pixel size, or a user-specified value. 
        
        After a 3D positional uncertainty sphere has been determined for each point in either dataframe, 
        the two orientations are compared, and points which overlap in 3D space from both orientations
        are saved as overlapped indexes.
        
        Attributes:
        z_range: int
            The radius of the Z undertainty for the 3D sphere for the XY data, which does not naturally
            include Z uncertainty values. 
            The recommended values are either the Z step [nm] used between successive images, or the 
            width of the light sheet, which designates the in-depth resolution of the system.
        xy_range: int
            The radius of the Y uncertainty for the 3D sphere for the XZ data, which does not naturally
            include Y uncertainty values.
            The recommended value is the in-plane pixel size [nm], which designates the Y step [nm]
            when the in-plane data is subjected to a pixelwise transformation without interpolation.
            
        Return:
            None. Will modify the data established in place.
        """
        range_x_xz = []
        range_x_xy = []
        range_y_xz = []
        range_y_xy = []
        range_z_xz = []
        range_z_xy = []
        for i in range(0, len(self.data.dfxz["X_XZ"])):
            x_left_xz = (self.data.dfxz["X_XZ"][i]-self.data.dfxz["U_X"][i])
            x_right_xz = (self.data.dfxz["X_XZ"][i]+self.data.dfxz["U_X"][i])
            y_left_xz = (self.data.dfxz["Y_XZ"][i]-xy_range)
            y_right_xz = (self.data.dfxz["Y_XZ"][i]+xy_range)
            z_left_xz = (self.data.dfxz["Z_XZ"][i]-self.data.dfxz["U_Z"][i])
            z_right_xz = (self.data.dfxz["Z_XZ"][i]+self.data.dfxz["U_Z"][i])
            range_x_xz.append((x_left_xz, x_right_xz))
            range_y_xz.append((y_left_xz, y_right_xz))
            range_z_xz.append((z_left_xz, z_right_xz))
        RX_xz = pd.arrays.IntervalArray.from_tuples(range_x_xz)
        RY_xz = pd.arrays.IntervalArray.from_tuples(range_y_xz)
        RZ_xz = pd.arrays.IntervalArray.from_tuples(range_z_xz)
        new_columns_data = {"X Range XZ": RX_xz, 
                            "Y Range XZ": RY_xz, 
                            "Z Range XZ": RZ_xz}  
        new_columns_df = pd.DataFrame(new_columns_data)
        self.data.dfxz = pd.concat([self.data.dfxz, new_columns_df], axis=1)
        for j in range(0, len(self.data.dfxy["X_XY"])):
            x_left_xy = (self.data.dfxy["X_XY"][j]-self.data.dfxy["U_XY"][j])
            x_right_xy = (self.data.dfxy["X_XY"][j]+self.data.dfxy["U_XY"][j])
            y_left_xy = (self.data.dfxy["Y_XY"][j]-self.data.dfxy["U_XY"][j])
            y_right_xy = (self.data.dfxy["Y_XY"][j]+self.data.dfxy["U_XY"][j])
            z_left_xy = (self.data.dfxy["Z_XY"][j]-z_range)
            z_right_xy = (self.data.dfxy["Z_XY"][j]+z_range)
            range_x_xy.append((x_left_xy, x_right_xy))
            range_y_xy.append((y_left_xy, y_right_xy))
            range_z_xy.append((z_left_xy, z_right_xy))
        RX_xy = pd.arrays.IntervalArray.from_tuples(range_x_xy)
        RY_xy = pd.arrays.IntervalArray.from_tuples(range_y_xy)
        RZ_xy = pd.arrays.IntervalArray.from_tuples(range_z_xy)
        new_columns_datay = {"X Range XY": RX_xy, 
                             "Y Range XY": RY_xy, 
                             "Z Range XY": RZ_xy}  
        new_columns_dfy = pd.DataFrame(new_columns_datay)
        self.data.dfxy = pd.concat([self.data.dfxy, new_columns_dfy], axis=1)
        all_indexes_XY = []
        all_indexes_XZ = []
        for k in range(0, len(self.data.dfxz["X_XZ"])):
            X_overlap = np.where(RX_xy.overlaps(RX_xz[k]))
            if X_overlap[0].size != 0:
                Y_overlap = np.where(RY_xy.overlaps(RY_xz[k]))
                xy = np.intersect1d(X_overlap, Y_overlap)
                if xy.size != 0:
                    Z_overlap = np.where(RZ_xy.overlaps(RZ_xz[k]))
                    xyz = np.intersect1d(xy, Z_overlap)
                    if xyz.size != 0:
                        value = [k]*len(xyz)
                        all_indexes_XY.append(xyz.tolist())
                        all_indexes_XZ.append(value)
        self.XY_indexes = sum(all_indexes_XY, []) 
        self.XZ_indexes = sum(all_indexes_XZ, [])
        return self
    
    def values(self):
        """
        A technique to determine the positional information using the indexes of overlap
        determined within the method.
        
        Attributes:
            None.
        Return:
            None. Will modify the data established in place. 
        """
        df = pd.DataFrame(columns=["X_XY", "Y_XY", "Z_XY", 
                                   "U_XY", "I_XY", "S_XY", 
                                   "X_XZ", "Y_XZ", "Z_XZ", 
                                   "U_X", "U_Z", "I_XZ", "S_X", "S_Z"])
        for i in range(0, len(self.XY_indexes)):
            X_XY = self.data.dfxy["X_XY"][self.XY_indexes[i]]
            Y_XY = self.data.dfxy["Y_XY"][self.XY_indexes[i]]
            Z_XY = self.data.dfxy["Z_XY"][self.XY_indexes[i]]
            U_XY = self.data.dfxy["U_XY"][self.XY_indexes[i]]
            I_XY = self.data.dfxy["I_XY"][self.XY_indexes[i]]
            S_XY = self.data.dfxy["S_XY"][self.XY_indexes[i]]
            X_XZ = self.data.dfxz["X_XZ"][self.XZ_indexes[i]]
            Y_XZ = self.data.dfxz["Y_XZ"][self.XZ_indexes[i]]
            Z_XZ = self.data.dfxz["Z_XZ"][self.XZ_indexes[i]]
            U_X = self.data.dfxz["U_X"][self.XZ_indexes[i]]
            U_Z = self.data.dfxz["U_Z"][self.XZ_indexes[i]]
            I_XZ = self.data.dfxz["I_XZ"][self.XZ_indexes[i]]
            S_X = self.data.dfxz["S_X"][self.XZ_indexes[i]]
            S_Z = self.data.dfxz["S_Z"][self.XZ_indexes[i]]
            df.loc[i] = [X_XY]+[Y_XY]+[Z_XY]+[U_XY]+[I_XY]+[S_XY]+[X_XZ]+[Y_XZ]+[Z_XZ]+[U_X]+[U_Z]+[I_XZ]+[S_X]+[S_Z]
        self.df=df
        self.dfxy = pd.concat([df["X_XY"], df["Y_XY"], df["Z_XY"], df["U_XY"], df["I_XY"], df["S_XY"]], 
                            keys=["X_XY", "Y_XY", "Z_XY", "U_XY", "I_XY", "S_XY"], axis=1)
        self.dfxz = pd.concat([df["X_XZ"], df["Y_XZ"], df["Z_XZ"], df["U_X"], df["U_Z"], df["I_XY"], df["S_X"], df["S_Z"]], 
                            keys=["X_XZ", "Y_XZ", "Z_XZ", "U_X", "U_Z", "I_XZ", "S_X", "S_Z"], axis=1)
        return self

    def download_dataframe(self):
        """
        A technique to download the data prepared by the Overlap.py method as a CSV file named
        "Overlap.py_dfXY_dfXZ.csv".
        
        Attributes:
            None.
        Return:
            None. Will download the dataframe as a CSV file.
        """
        download_df = pd.concat([self.dfxy, self.dfxz], axis=1, sort=False)
        download_df = download_df.rename(columns={'X_XY': 'x_xy [nm]', 
                                                  'Y_XY': 'y_xy [nm]', 
                                                  'Z_XY': 'z_xy [nm]', 
                                                  'U_XY': 'uncertainty_xy [nm]', 
                                                  'I_XY': 'intensity_xy [photon]', 
                                                  'S_XY': 'sigma_xy [nm]', 
                                                  'X_XZ': 'x_xz [nm]', 
                                                  'Y_XZ': 'y_xz [nm]', 
                                                  'Z_XZ': 'z_xz [nm]', 
                                                  'U_X': 'uncertainty_x [nm]',
                                                  'U_Z': 'uncertainty_z [nm]',
                                                  'I_XZ': 'intensity_xz [photon]',
                                                  'S_X': 'sigma_x [nm]', 
                                                  'S_Z': 'sigma_z [nm]'})
        download_df.to_csv("Overlap.py_dfXY_dfXZ.csv", index=False, encoding='utf-8')
        return self