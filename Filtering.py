class filtering(object):
    """
    Sierra Dean
    RIKEN SPring-8 Center
    March 4 2024
    """
    def __init__(self, data):
        """
        A technique to filter all of the points overlapping in 3D space.
        
        Attributes:
        data: 
            The data previously developed and contained within the overlap class.
 
        Return:
            None. Will modify the data established in place.
        """
        self.data = data
        return 
    
    def merge(self):
        """
        A technique to determine all of the indexes for points which overlap with multiple other points.
 
        Return:
            None. Will modify the data established in place.
        """
        unique_XY = np.unique(self.data.df["X_XY"])
        unique_XZ = np.unique(self.data.df["X_XZ"])
        XY_indexes = []
        XZ_indexes = []
        for i in unique_XZ:
            XZ_indexes.append(np.where(self.data.df["X_XZ"]==i)[0].tolist())
        for j in unique_XY:
            XY_indexes.append(np.where(self.data.df["X_XY"]==j)[0].tolist())
        all_items = np.arange(0, len(self.data.df["X_XZ"]))
        all_indexes = []
        for k in all_items:
            list_all = np.array([k])
            for l in list_all:
                list_z = []
                list_y = []
                [list_z.append(sub_list) for sub_list in XZ_indexes if l in sub_list]
                [list_y.append(sub_list) for sub_list in XY_indexes if l in sub_list]
                list_all = np.unique(np.append(list_all, list_z))
                list_all = np.unique(np.append(list_all, list_y))
            all_indexes.append(list_all.tolist())
        self.all_indexes = all_indexes
        merged_indexes = []
        for m in all_indexes:
            m = set(m)
            for n in merged_indexes:
                if m & n:
                    n.update(m)
                    break
            else:
                merged_indexes.append(m)
        merged_indexes = [list(o) for o in merged_indexes]
        self.merged_indexes = merged_indexes
        return self
    
    def selection(self, selection_type="uncertainty"):
        """
        A technique to determine all of the indexes for points which overlap with multiple other points.
        
        Attributes:
        selection_type: str "uncertainty", "Uncertainty", "intensity", "Intensity"
            The method of filtering. If uncertainty is selected, for all overlapped points, 
            only the lowest positional uncertainty point will be kept. If intensity is selected, 
            only the lowest intensity point will be kept.
 
        Return:
            None. Will modify the data established in place.
        """
        if (selection_type=="uncertainty") or (selection_type=="Uncertainty"):
            point_selection = []
            for i in range(0, len(self.merged_indexes)):
                u_xy = []
                u_z = []
                u_x = []
                sumall = []
                for j in self.merged_indexes[i]:
                    u_xy.append(self.data.df["U_XY"][j])
                    u_z.append(self.data.df["U_Z"][j])
                    u_x.append(self.data.df["U_X"][j])
                [sumall.append(sum(s)) for s in zip(*[u_xy, u_z, u_x])]
                loc = np.where(sumall==min(sumall))[0]
                point_selection.append(self.merged_indexes[i][loc[0]])
            self.point_indexes = point_selection
        elif (selection_type=="intensity") or (selection_type=="Intensity"):
            point_selection = []
            for i in range(0, len(self.merged_indexes)):
                i_xy = []
                i_xz = []
                sumall = []
                for j in self.merged_indexes[i]:
                    i_xy.append(self.data.df["I_XY"][j])
                    i_xz.append(self.data.df["I_XZ"][j])
                [sumall.append(sum(s)) for s in zip(*[i_xy, i_xz])]
                loc = np.where(sumall==min(sumall))[0]
                point_selection.append(self.merged_indexes[i][loc[0]])
            self.point_indexes = point_selection
        return self
    
    def points(self):
        """
        A technique to obtain the point information from the indexes obtained from filtering.

        Return:
            None. Will modify the data established in place.
        """
        df = self.data.df
        for i in range(0, len(df)):
            if i not in (self.point_indexes):
                df = df.drop(i, axis=0)
        self.df = df
        points = pd.concat([df["X_XY"], df["Y_XY"], df["Z_XZ"], 
                           df["U_XY"], df["U_Z"], 
                           df["S_XY"], df["S_Z"], 
                            df["I_XY"], df["I_XZ"]], 
                           keys = ["X [nm]", "Y [nm]", "Z [nm]", 
                                   "Uncertainty XY [nm]", "Uncertainty Z [nm]",
                                   "Sigma XY [nm]", "Sigma Z [nm]",
                                   "Intensity XY [Photons]", "Intensity XZ [Photons]"], axis=1).reset_index()
        self.points = points
        return self