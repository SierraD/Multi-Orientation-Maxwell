class surface(object):
    """
    This file is part of the 3D STORM software
    
    File author(s): Sierra Dean <ccnd@live.com>
    
    Distributed under the GPLv3 Licence.
    See accompanying file LICENSE.txt or copy at
        http://www.gnu.org/licenses/gpl-3.0.html
        
    source: https://github.com/SierraD/3DSTORM
    
    Last Updated: June 18 2024
    """
    def __init__(self, data):
        """
        A technique to estimate the radius and center position of a sphere using a linear 
        least-square sphere-fitting algorithm, namely a Summation Least-Squares (SLS) fitting.
        
        The three-dimensional data returned by the Filtering.py class (i.e. X [nm], 
        Uncertainty XY [nm], etc.) is first translated to a new coordinate space throgh domain 
        shifting (x, y, z)->(u, v, w), with the least-squares minimising function then linearised,
        and partially differentiated for minimization before translating back to the original 
        coordinate system.
        
        Attributes: 
        data: 
            The data previously developed and contained within the filtering class.
            
        Return: 
            None. Will modify the data established in place. 
        """
        self.data = data
        return
    
    def evaluation(self):
        """
        A technique to determine the radius and center position of the data returned 
        by the Filtering class using a Summation Least-Squares sphere fitting.
        
        Attributes: 
            None.
            
        Return: 
            None. Will modify the data established in place. 
        """   
        N = len(self.data.points["X [nm]"])
        
        u_i = self.data.points["X [nm]"] - sum(self.data.points["X [nm]"])/N
        v_i = self.data.points["Y [nm]"] - sum(self.data.points["Y [nm]"])/N
        w_i = self.data.points["Z [nm]"] - sum(self.data.points["Z [nm]"])/N
        
        s_uu = sum(u_i * u_i)
        s_uv = sum(u_i * v_i)
        s_uw = sum(u_i * w_i)
        s_vv = sum(v_i * v_i)
        s_vw = sum(v_i * w_i)
        s_ww = sum(w_i * w_i)
        
        s_uuu = sum(u_i * u_i * u_i)
        s_uvv = sum(u_i * v_i * v_i)
        s_uww = sum(u_i * w_i * w_i)
        s_uuv = sum(u_i * u_i * v_i)
        s_vvv = sum(v_i * v_i * v_i)
        s_vww = sum(v_i * w_i * w_i)
        s_uuw = sum(u_i * u_i * w_i)
        s_vvw = sum(v_i * v_i * w_i)
        s_www = sum(w_i * w_i * w_i)
        
        A = [[s_uu, s_uv, s_uw], [s_uv, s_vv, s_vw], [s_uw, s_vw, s_ww]]
        B = [[s_uuu+s_uvv+s_uww], [s_uuv+s_vvv+s_vww], [s_uuw+s_vvw+s_www]]
        x = 0.5*numpy.matmul(numpy.linalg.inv(A), B)
        
        self.X_cent = [x[0] + sum(self.data.points["X [nm]"])/N][0][0]
        self.Y_cent = [x[1] + sum(self.data.points["Y [nm]"])/N][0][0]
        self.Z_cent = [x[2] + sum(self.data.points["Z [nm]"])/N][0][0]
        
        self.radius = numpy.sqrt(x[0]**2 + x[1]**2 + x[2]**2 + ((s_uu+s_vv+s_ww)/N))
        self.radius_i = numpy.sqrt((self.data.points["X [nm]"] - self.X_cent)**2 + (self.data.points["Y [nm]"] - self.Y_cent)**2 + (self.data.points["Z [nm]"] - self.Z_cent)**2)
        self.error_i = self.radius_i - self.radius
        return self
