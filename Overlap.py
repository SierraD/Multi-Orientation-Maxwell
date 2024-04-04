class overlap(object):
    """
    Sierra Dean
    RIKEN SPring-8 Center
    March 4 2024
    """
    
    def __init__(self, data):
        """
        A technique to use two different orientations of ThunderSTORM results to determine the overlapping points in 3D space.
        
        Attributes:
        data: 
            The data previously developed and contained within the preparations class.
            
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
        A technique to determine the indexes of the dataframe where points overlap. This version, computationally simpler
        than its counterpart, defines a 3D sphere around each XY orientation point comprised of that point's uncertainty
        values, and searches to see if any XZ point exists within that 3D range. Although equally accurate, this method has 
        a tendency to dismiss particle overlaps since XZ is regarded as a ground truth data, and the XZ uncertainties are 
        not taken into account. Therefore, if two particles only overlap due to the positional uncertainty from both XY
        and XZ orientations, they will not be included.
        
        Attributes:
        z_range: int
            The search radius for the overlap in the Z direction, in nanometers. 
            Generally selected to be either the step size between images or the width of the lightsheet, if known.
            Represents the intial positional uncertainty of the Z dimension obtained from the XY orientation.
        xy_range: int
            The search radius for the overlap in the XY directions, in nanometers.
            Generally selected to be the image pixel size in the XY orientation.
            Represents the initial positional uncertainty of the Y dimension,
            obtained from rotating the XY data to the XZ orientation through a pixelwise transformation without interpolation.
            
        Return:
            None. Will modify the data established in place.
        """
        all_indexes_XY = []
        all_indexes_XZ = []
        for i in range(0, len(self.data.dfxz["X_XZ"])):
            x_left = np.where(self.data.dfxz["X_XZ"][i]-self.data.dfxz["U_X"][i]<self.data.dfxy["X_XY"])
            x_right = np.where(self.data.dfxy["X_XY"]<self.data.dfxz["X_XZ"][i]+self.data.dfxz["U_X"][i])
            x = np.intersect1d(x_left, x_right)
            y_left = np.where(self.data.dfxz["Y_XZ"][i]-xy_range<self.data.dfxy["Y_XY"])
            y_right = np.where(self.data.dfxy["Y_XY"]<self.data.dfxz["Y_XZ"][i]+xy_range)
            y = np.intersect1d(y_left, y_right)
            xy = np.intersect1d(x, y)
            if xy.size != 0:
                z_left = np.where(self.data.dfxz["Z_XZ"][i]-z_range<self.data.dfxy["Z_XY"])
                z_right = np.where(self.data.dfxy["Z_XY"]<self.data.dfxz["Z_XZ"][i]+z_range)
                z = np.intersect1d(z_left, z_right)
                all_pts = np.intersect1d(xy, z)
                if all_pts.size !=0:
                    value = [i]*len(all_pts)
                    all_indexes_XY.append(all_pts.tolist())
                    all_indexes_XZ.append(value)
        self.XY_indexes = sum(all_indexes_XY, []) 
        self.XZ_indexes = sum(all_indexes_XZ, [])
        return self
    
    def indexes_computational(self, z_range, xy_range):
        """
        A technique to determine the indexes of the dataframe where points overlap. This version, computationally heavy
        compared to its counterpart, defines 3D spheres around all points comprised of that point's uncertainty values, 
        and searches to see if any other 3D sphere exist within the defined range. Due to the computational cost of this
        method, it is recommended that the simplified version be run first, and this method only used if it is imperative
        to obtain overlap from all points.
        
        Attributes:
        z_range: int
            The search radius for the overlap in the Z direction, in nanometers. 
            Generally selected to be either the step size between images or the width of the lightsheet, if known.
            Represents the intial positional uncertainty of the Z dimension obtained from the XY orientation.
        xy_range: int
            The search radius for the overlap in the XY directions, in nanometers.
            Generally selected to be the image pixel size in the XY orientation.
            Represents the initial positional uncertainty of the Y dimension,
            obtained from rotating the XY data to the XZ orientation through a pixelwise transformation without interpolation.
            
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
        A technique to determine the indexes of the dataframe where points overlap.
        
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