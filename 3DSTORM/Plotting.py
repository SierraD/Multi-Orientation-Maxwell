class plotting(object):
    """
    Sierra Dean
    RIKEN SPring-8 Center
    March 4 2024
    """
    def __init__(self, data):
        """
        A technique to visualize the 2D and 3D results obtained from the method.
        
        Attributes:
        data: 
            The data previously developed and contained within the various classes.
 
        Return:
            A plot displaying the values as specified.
        """
        self.data = data
        return
    
    def tricolumn(self, dimensions=2):
        """
        A technique to visualize the 3D results as three columns of two dimensional data, 
        inlcuding the error bars for the localization.
        
        Attributes:
        dimensions: int 2 or 3 
            The dimensions of the data, either 2D for the XY/XZ data or 3D for the finalized data.
 
        Return:
            A plot displaying the localizations as three columns of 2D data.
        """
        fig = make_subplots(rows=1, cols=3, subplot_titles=("XY Orientation", "XZ Orientation", "YZ Orientation"))
        if dimensions == 2:
            fig.add_trace(go.Scatter(x=self.data.dfxy["X_XY"], y=self.data.dfxy["Y_XY"], name="XY",
                                     error_x=dict(type='data', array=self.data.dfxy["U_XY"], visible=True, width=1, color="gray"), 
                                     error_y=dict(type='data', array=self.data.dfxy["U_XY"], visible=True, width=1, color="gray"),   
                                     marker=dict(color="Red", size=2), mode="markers"), row=1, col=1)
            fig.add_trace(go.Scatter(x=self.data.dfxz["X_XZ"], y=self.data.dfxz["Y_XZ"], name="XZ",
                                     error_x=dict(type='data', array=self.data.dfxz["U_X"], visible=True, width=1, color="gray"), 
                                 marker=dict(color="Blue", size=2), mode="markers"), row=1, col=1)
            fig.add_trace(go.Scatter(x=self.data.dfxy["X_XY"], y=self.data.dfxy["Z_XY"], 
                                     error_x=dict(type='data', array=self.data.dfxy["U_XY"], visible=True, width=1, color="gray"), 
                                     marker=dict(color="Red", size=2), mode="markers", showlegend=False), row=1, col=2)
            fig.add_trace(go.Scatter(x=self.data.dfxz["X_XZ"], y=self.data.dfxz["Z_XZ"], 
                                     error_x=dict(type='data', array=self.data.dfxz["U_X"], visible=True, width=1, color="gray"), 
                                     error_y=dict(type='data', array=self.data.dfxz["U_Z"], visible=True, width=1, color="gray"),    
                                     marker=dict(color="Blue", size=2), mode="markers", showlegend=False), row=1, col=2)
            fig.add_trace(go.Scatter(x=self.data.dfxy["Y_XY"], y=self.data.dfxy["Z_XY"], 
                                     error_x=dict(type='data', array=self.data.dfxy["U_XY"], visible=True, width=1, color="gray"), 
                                     marker=dict(color="Red", size=2), mode="markers", showlegend=False), row=1, col=3)
            fig.add_trace(go.Scatter(x=self.data.dfxz["Y_XZ"], y=self.data.dfxz["Z_XZ"], 
                                     error_y=dict(type='data', array=self.data.dfxz["U_Z"], visible=True, width=1, color="gray"),    
                                     marker=dict(color="Blue", size=2), mode="markers", showlegend=False), row=1, col=3)
        elif dimensions == 3:
            fig.add_trace(go.Scatter(x=self.data.points["X [nm]"], y=self.data.points["Y [nm]"],
                                     error_x=dict(type='data', array=self.data.points["Uncertainty XY [nm]"], visible=True, width=1, color="gray"), 
                                     error_y=dict(type='data', array=self.data.points["Uncertainty XY [nm]"], visible=True, width=1, color="gray"),   
                                     marker=dict(color="Black", size=2), mode="markers", showlegend=False), row=1, col=1)
            fig.add_trace(go.Scatter(x=self.data.points["X [nm]"], y=self.data.points["Z [nm]"],
                                     error_x=dict(type='data', array=self.data.points["Uncertainty XY [nm]"], visible=True, width=1, color="gray"), 
                                     error_y=dict(type='data', array=self.data.points["Uncertainty Z [nm]"], visible=True, width=1, color="gray"),   
                                     marker=dict(color="Black", size=2), mode="markers", showlegend=False), row=1, col=2)
            fig.add_trace(go.Scatter(x=self.data.points["Y [nm]"], y=self.data.points["Z [nm]"],
                                     error_x=dict(type='data', array=self.data.points["Uncertainty XY [nm]"], visible=True, width=1, color="gray"), 
                                     error_y=dict(type='data', array=self.data.points["Uncertainty Z [nm]"], visible=True, width=1, color="gray"),   
                                     marker=dict(color="Black", size=2), mode="markers", showlegend=False), row=1, col=3)
        fig.update_xaxes(title_text="X [nm]", row=1, col=1, title_standoff = 0)
        fig.update_yaxes(title_text="Y [nm]", row=1, col=1, title_standoff = 0)
        fig.update_xaxes(title_text="X [nm]", row=1, col=2, title_standoff = 0)
        fig.update_yaxes(title_text="Z [nm]", row=1, col=2, title_standoff = 0)
        fig.update_xaxes(title_text="Y [nm]", row=1, col=3, title_standoff = 0)
        fig.update_yaxes(title_text="Z [nm]", row=1, col=3, title_standoff = 0)
        fig.show()
        return self
    
    def tricolumn_sigma(self, dimensions=2):
        """
        A technique to visualize the 3D results as three columns of two dimensional data, 
        with the size of the markers represented as the sigma value obtained from ThunderSTORM.
        
        Attributes:
        dimensions: int 2 or 3 
            The dimensions of the data, either 2D for the XY/XZ data or 3D for the finalized data.
 
        Return:
            A plot displaying the localizations as three columns of 2D data visualized as the ThunderSTORM sigma value.
        """
        if dimensions == 2:
            fig = make_subplots(rows=1, cols=2, subplot_titles=("XY Orientation", "XZ Orientation"))
            left_x_xy = self.data.dfxy["X_XY"] - self.data.dfxy["S_XY"]
            right_x_xy = self.data.dfxy["X_XY"] + self.data.dfxy["S_XY"]
            left_y = self.data.dfxy["Y_XY"] - self.data.dfxy["S_XY"]
            right_y = self.data.dfxy["Y_XY"] + self.data.dfxy["S_XY"]
            left_x_xz = self.data.dfxz["X_XZ"] - self.data.dfxz["S_X"]
            right_x_xz = self.data.dfxz["X_XZ"] + self.data.dfxz["S_X"]
            left_z = self.data.dfxz["Z_XZ"] - self.data.dfxz["S_Z"]
            right_z = self.data.dfxz["Z_XZ"] + self.data.dfxz["S_Z"]
            for i in range(0, len(left_x_xy)): 
                fig.add_shape(type="circle", x0=left_x_xy[i], x1=right_x_xy[i], y0=left_y[i], y1=right_y[i], 
                              fillcolor="red", line_color="red", opacity=0.2, row=1, col=1)
            for j in range(0, len(left_x_xz)): 
                fig.add_shape(type="circle", x0=left_x_xz[j], x1=right_x_xz[j], y0=left_z[j], y1=right_z[j], 
                              fillcolor="Blue", line_color="Blue", opacity=0.2, row=1, col=2)
            fig.update_xaxes(range=[min(left_x_xy)+((min(right_x_xy)-max(left_x_xy))/4),
                                max(right_x_xy)-((min(right_x_xy)-max(left_x_xy))/4)], row=1, col=1)
            fig.update_yaxes(range=[min(left_y)+((min(right_y)-max(left_y))/4),
                                max(right_y)-((min(right_y)-max(left_y))/4)], row=1, col=1)
            fig.update_xaxes(range=[min(left_x_xz)+((min(right_x_xz)-max(left_x_xz))/4),
                                max(right_x_xz)-((min(right_x_xz)-max(left_x_xz))/4)], row=1, col=2)
            fig.update_yaxes(range=[min(left_z)+((min(right_z)-max(left_z))/4),
                                max(right_z)-((min(right_z)-max(left_z))/4)], row=1, col=2)
            fig.show()
        elif dimensions == 3:
            fig = make_subplots(rows=1, cols=3, subplot_titles=("XY Orientation", "XZ Orientation", "YZ Orientation"))
            left_x = self.data.points["X [nm]"] - self.data.points["Sigma XY [nm]"]
            right_x = self.data.points["X [nm]"] + self.data.points["Sigma XY [nm]"]
            left_y = self.data.points["Y [nm]"] - self.data.points["Sigma XY [nm]"]
            right_y = self.data.points["Y [nm]"] + self.data.points["Sigma XY [nm]"]
            left_z = self.data.points["Z [nm]"] - self.data.points["Sigma Z [nm]"]
            right_z = self.data.points["Z [nm]"] + self.data.points["Sigma Z [nm]"]
            for i in range(0, len(left_x)): 
                fig.add_shape(type="circle", x0=left_x[i], x1=right_x[i], y0=left_y[i], y1=right_y[i], 
                              fillcolor="Green", line_color="Green", opacity=0.2, row=1, col=1)
                fig.add_shape(type="circle", x0=left_x[i], x1=right_x[i], y0=left_z[i], y1=right_z[i], 
                              fillcolor="Green", line_color="Green", opacity=0.2, row=1, col=2)
                fig.add_shape(type="circle", x0=left_y[i], x1=right_y[i], y0=left_z[i], y1=right_z[i], 
                              fillcolor="Green", line_color="Green", opacity=0.2, row=1, col=3)
            fig.update_xaxes(range=[min(left_x)+((min(right_x)-max(left_x))/4),
                                max(right_x)-((min(right_x)-max(left_x))/4)], row=1, col=1)
            fig.update_yaxes(range=[min(left_y)+((min(right_y)-max(left_y))/4),
                                max(right_y)-((min(right_y)-max(left_y))/4)], row=1, col=1)
            fig.update_xaxes(range=[min(left_x)+((min(right_x)-max(left_x))/4),
                                max(right_x)-((min(right_x)-max(left_x))/4)], row=1, col=2)
            fig.update_yaxes(range=[min(left_z)+((min(right_z)-max(left_z))/4),
                                max(right_z)-((min(right_z)-max(left_z))/4)], row=1, col=2)
            fig.update_xaxes(range=[min(left_y)+((min(right_y)-max(left_y))/4),
                                max(right_y)-((min(right_y)-max(left_y))/4)], row=1, col=3)
            fig.update_yaxes(range=[min(left_z)+((min(right_z)-max(left_z))/4),
                                max(right_z)-((min(right_z)-max(left_z))/4)], row=1, col=3)
            fig.show()
        return
    
    def three_dimensional(self, dimensions=2):
        """
        A technique to visualize the 3D results of the localizations.
        
        Attributes:
        dimensions: int 2 or 3 
            The dimensions of the data, either 2D for the XY/XZ data or 3D for the finalized data.
 
        Return:
            A 3D visualization of the localizations.
        """
        fig = go.Figure()
        if dimensions == 2:
            fig.add_trace(go.Scatter3d(x=self.data.dfxy["X_XY"], y=self.data.dfxy["Y_XY"], z=self.data.dfxy["Z_XY"],
                                       marker=dict(color="#FF2D00", size=2), mode="markers"))
            fig.add_trace(go.Scatter3d(x=self.data.dfxz["X_XZ"], y=self.data.dfxz["Y_XZ"], z=self.data.dfxz["Z_XZ"],
                                       marker=dict(color="#001BFF", size=2), mode="markers"))
        elif dimensions == 3:
            fig.add_trace(go.Scatter3d(x=self.data.points["X [nm]"], y=self.data.points["Y [nm]"], z=self.data.points["Z [nm]"],
                                       marker=dict(color="#000000", size=2), mode="markers"))
        fig.show()
        return
    
    def three_dimensional_sigma(self, dimensions=3):
        """
        A technique to visualize the 3D results of the localizations, 
        with the size represented as the sigma value obtained from ThunderSTORM.
        
        Attributes:
        dimensions: int 2
            The dimensions of the data, should be 3D for the finalized data.
 
        Return:
            A 3D visualization of the localizations with the size represented as the 
            sigma value obtained from ThunderSTORM.
        """
        if dimensions != 3:
            raise ValueError("The PSF is only three dimensional after the data has been converted to three dimensions.")
        u, v = np.mgrid[0:2*np.pi:20j, 0:np.pi:10j]
        fig = make_subplots(rows=1, cols=1, specs=[[{'is_3d': True}]])
        for i in range(0,len(self.data.points["X [nm]"])):
            x = self.data.points["Sigma XY [nm]"][i]*np.cos(u)*np.sin(v)+self.data.points["X [nm]"][i]
            y = self.data.points["Sigma XY [nm]"][i]*np.sin(u)*np.sin(v)+self.data.points["Y [nm]"][i]
            z = self.data.points["Sigma Z [nm]"][i]*np.cos(v)+self.data.points["Z [nm]"][i]
            fig.add_trace(go.Surface(x=x, y=y, z=z, opacity=0.5), 1, 1)
        fig.update_traces(showscale=False)
        fig.update_layout(scene = dict(xaxis_title="X [nm]", yaxis_title="Y [nm]", zaxis_title="Z [nm]"))
        fig.show()
        return
