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
        A technique to determine the indexes of the dataframe where points overlap.
        
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
