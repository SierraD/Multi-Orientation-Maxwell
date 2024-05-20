class precision(object):
    """
    This file is part of the 3D STORM software
    
    File author(s): Sierra Dean <ccnd@live.com>
    
    Distributed under the GPLv3 Licence.
    See accompanying file LICENSE.txt or copy at
        http://www.gnu.org/licenses/gpl-3.0.html
        
    source: https://github.com/SierraD/3DSTORM
    
    Last Updated: May 20 2024
    """
    def __init__(self, data, data_type="Uncertainty XY [nm]",):
        """
        A technique to automatically obtain the precision of the method. 
        
        The data returned by ThunderSTORM (i.e. Uncertainty, Sigma) are standard deviation 
        values obtained by fitting the fluorescent spot many times, and therefore the sigma
        value obtained from this deviation can be interpreted as the method precision. 
        
        To calculate the method precision, the histogram of the data is automatically 
        compared with Probability Distribution Functions (PDF) from known distributions 
        (i.e. Normal, Chi2, Exponential, etc.), then the best fit is selected by minimizing 
        the Sum of the Square Error between the histogram and the PDF.
        
        
        Attributes:
        data: 
            The data previously developed and contained within the filtering class.
 
        Return:
            None. Will modify the data established in place.
        """
        self.data = data
        self.data_type = data_type
        self.points = self.data.points
        self.height = data.points[data_type].values
        self.best_fit = "norm"
        return
    
    def get_precision(self, distribution=get_common_distributions()):
        """
        A technique to generate a histogram from specified data, and compare with PDF fittings from 
        known distributions. 
        
        A plot of the data is automatically generated, and up to five PDF distributions, selected 
        by minizming the Sum of the Square Error, are plotted in conjunction.
                
        Attributes:
        data: str "Uncertainty XY [nm]", "Uncertainty Z [nm]", "Sigma XY [nm]", "Sigma Z [nm]", etc.
            The data to be used for the histogram and the fitting.
        distribution: str "norm", "lognorm", get_common_distributions(), get_distributions(), etc.
            The type of distribution to be used for the fitting.
            If none specified, all common distributions (10) will be compared.
            
        Return:
            None. Will modify the data established in place.
        """
        fit = Fitter(self.height, distributions=distribution)
        fit.fit()
        print(fit.summary())
        self.best_fit = list(fit.get_best().keys())[0]
        print("\nBest Fit Selected: ", self.best_fit)
        return self
    
    def fitted_precision(self, kde=False, bins=30, distribution=False):
        """
        A technique to automatically obtain the precision of data using either the determined or 
        specified best fit from known distributions. 
        
        A plot of the data is automatically generated, with the specified PDF plotted in 
        conjunction, in addition to displaying useful values related to method precision.
        
        Attributes:
        kde: bool 
            The decision to also display the Kernel Density Estimation (KDE) of the histogram.
        bins: int
            The number of bins to be used for the final histogram.
        distribution: False or str "norm", "lognorm", etc.
            Eiter False, which will automatically use the best fit determined by the 
            'get_precision' method, or a string specifying the known distribution to 
            perform a PDF comparison with.

        Return:
            A plot displaying the histogram and important variables.
        """
        if distribution == False:
            distribution = self.best_fit
        else: 
            distribution = distribution
        
        fit = Fitter(self.height, distributions=distribution)
        fit.fit()
        fit.summary()
        xmin, xmax = plt.xlim()
        self.Xaxis = np.linspace(xmin, xmax, 1000)
        
        params = list(fit.fitted_param[distribution]) 
        best_fitted = eval("scipy.stats." + distribution)
        fitting = best_fitted.pdf(self.Xaxis, *params)
        
        fig0, ax0 = plt.subplots(figsize=(10,5))
        ax0.set_title("Precision Analysis")
        ax0.set_ylabel("Density")
        ax0.set_xlabel(self.data_type)
        sns.histplot(x=self.points[self.data_type], ax=ax0, stat="density", kde=kde, bins=bins, linewidth=0.5, 
                     label=("N="+str(len(self.points[self.data_type]))), 
                     line_kws={'label': "Kernel Density Estimation (KDE) Fitting"})
        ax0.plot(self.Xaxis, fitting, color="red", linestyle="dashed", 
                 label="Probability Density Function (PDF) Fitting")
        a = round(self.Xaxis[np.where(fitting == fitting.max())][0], 4)
        b = round(fitting[np.where(fitting == fitting.max())][0], 4)
        ax0.scatter(a, b, color="red", label=("["+str(a)+", "+str(b)+"]"))
        ax0.legend()
        plt.show()
        return self