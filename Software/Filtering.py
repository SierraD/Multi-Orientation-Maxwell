class filtering(object):
    """
    This file is part of the Multi-Orientation MAXWELL software
    
    File author(s): Sierra Dean <ccnd@live.com>
    
    Distributed under the GPLv3 Licence.
    See accompanying file LICENSE.txt or copy at
        http://www.gnu.org/licenses/gpl-3.0.html
        
    source: https://github.com/SierraD/Multi-Orientation-MAXWELL
    
    Last Updated: Sept 26 2024
    """
    def __init__(self, data):
        """
        A technique to obtain precise localizations in 3D space, with data originating 
        from two-dimensional ThunderSTORM analysis performed on mulitple orientations. 
        
        This method uses the overlapped localizations determined by the Overlap.py method and 
        filters all overlapped points to ensure no duplications. The choice of filtering
        includes selecting the lowest uncertainty values, or the most similar intensity
        values. 
        
        The final dataframe obtained from this method is the resulting 3D localizations
        obtained from comparing the two ThunderSTORM orientations.
        
        Attributes:
        data: 
            The data previously developed and contained within the Overlap.py class.
 
        Return:
            None. Will modify the data established in place.
        """
        self.data = data
        return 
    
    def merge(self):
        """
        A technique to group indexes which contain points which overlap with multiple other
        points, which can then be used to filter out duplicates.
        
        Attributes:
            None. 
 
        Return:
            None. Will modify the data established in place.
        """
        unique_XY = numpy.unique(self.data.df["X_XY"])
        unique_XZ = numpy.unique(self.data.df["X_XZ"])
        XY_indexes = []
        XZ_indexes = []
        for i in unique_XZ:
            XZ_indexes.append(numpy.where(self.data.df["X_XZ"]==i)[0].tolist())
        for j in unique_XY:
            XY_indexes.append(numpy.where(self.data.df["X_XY"]==j)[0].tolist())
        all_items = numpy.arange(0, len(self.data.df["X_XZ"]))
        all_indexes = []
        for k in all_items:
            list_all = numpy.array([k])
            for l in list_all:
                list_z = []
                list_y = []
                [list_z.append(sub_list) for sub_list in XZ_indexes if l in sub_list]
                [list_y.append(sub_list) for sub_list in XY_indexes if l in sub_list]
                list_all = numpy.unique(numpy.append(list_all, list_z))
                list_all = numpy.unique(numpy.append(list_all, list_y))
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
        A technique to filter all of the indexes which contain overlaps from multiple 
        localizations to remove all non-unique localizations.
                
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
                loc = numpy.where(sumall==min(sumall))[0]
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
                loc = numpy.where(sumall==min(sumall))[0]
                point_selection.append(self.merged_indexes[i][loc[0]])
            self.point_indexes = point_selection
        return self
    
    def points(self):
        """
        A technique to obtain a dataframe of precise localizations in 3D space, 
        with data originating from two-dimensional ThunderSTORM analysis performed 
        in multiple orientations.
        
        Attributes:
            None.
            
        Return:
            None. Will modify the data established in place.
        """
        df = self.data.df
        for i in range(0, len(df)):
            if i not in (self.point_indexes):
                df = df.drop(i, axis=0)
        self.df = df
        points = pandas.concat([df["X_XY"], df["Y_XY"], df["Z_XZ"], 
                            df["U_XY"], df["U_Z"], 
                            df["S_XY"], df["S_Z"], 
                            df["I_XY"], df["I_XZ"], 
                            df["O_XY"], df["O_XZ"],
                            df["B_XY"], df["B_XZ"]], 
                           keys = ["X [nm]", "Y [nm]", "Z [nm]", 
                                   "Uncertainty XY [nm]", "Uncertainty Z [nm]",
                                   "Sigma XY [nm]", "Sigma Z [nm]",
                                   "Intensity XY [Photons]", "Intensity XZ [Photons]",
                                   "Offset XY [Photons]", "Offset XZ [Photons]",
                                   "Bkgstd XY [Photons]", "Bkgstd XZ [Photons]"], axis=1).reset_index(drop=True)
        self.points = points
        return self
    
    def download_dataframe(self, filename="3DSTORM"):
        """
        A technique to download the data prepared by the Filtering.py method as a 
        CSV file named "3DSTORM.csv".
        
        Attributes:
            None.
            
        Return:
            None. Will download the dataframe as a CSV file.
        """
        self.points.index.set_names('id', level=None, inplace=True)
        self.points.to_csv(filename+".csv", index=True, encoding='utf-8')
        return self
