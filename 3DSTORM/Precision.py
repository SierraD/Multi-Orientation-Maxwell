class final(object):
    """
    Sierra Dean
    RIKEN SPring-8 Center
    March 4 2024
    """
    def __init__(self, data):
        """
        A technique to obtain the 3D localization after filtering.
        
        Attributes:
        data: 
            The data previously developed and contained within the filtering class.
 
        Return:
            None. Will modify the data established in place.
        """
        self.data = data
        return
    
    def points(self):
        """
        A technique to obtain the point information from the indexes obtained from filtering.

        Return:
            None. Will modify the data established in place.
        """
        df = self.data.df
        for i in range(0, len(df)):
            if i not in (self.data.point_indexes):
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
    
    def get_distribution(self, data="Uncertainty XY [nm]", distribution=get_common_distributions(), method="sumsquare_error"):
        """
        A technique to provide a histogram for specified data and calculate a Kernel Density Estimation, 
        which is then fitted by common distributions to find the best fit for the data.
        
        Attributes:
        data: str "Uncertainty XY [nm]", "Uncertainty Z [nm]", "Sigma XY [nm]", "Sigma Z [nm]", etc
            The data to be used for the histogram and the fitting.
        distribution: str "norm", "lognorm", "gamma", "chi2", "exponpow", "rayleigh"
            The type of distribution to be used for the fitting.
            If none specified, all distributions will be compared for the best fit, 
            with the best fitting selected by the method.
            If selected, the variables will be saved to be used to calculate the precision from the fitting.
        method: str "sumsquare_error", etc
            The method of determining the best fit. Can be selected from the fit options.

        Return:
            None. Will modify the data established in place.
        """
        self.data = data
        self.distribution = distribution
        height = self.points[data].values
        fit = Fitter(height, distributions=distribution)
        fit.fit()
        print(fit.summary())
        print(fit.get_best(method=method))
        xmin, xmax = plt.xlim()
        self.Xaxis = np.linspace(xmin, xmax, 1000)
        if distribution == "norm":
            self.distribution = "norm"
            self.loc = fit.fitted_param[distribution][0]
            self.scale = fit.fitted_param[distribution][1]
        elif distribution == "lognorm":
            self.distribution = "lognorm"
            self.s = fit.fitted_param[distribution][0]
            self.loc = fit.fitted_param[distribution][1]
            self.scale = fit.fitted_param[distribution][2]
        elif distribution == "gamma":
            self.distribution = "gamma"
            self.a = fit.fitted_param[distribution][0]
            self.loc = fit.fitted_param[distribution][1]
            self.scale = fit.fitted_param[distribution][2]
        elif distribution == "chi2":
            self.distribution = "chi2"
            self.df = fit.fitted_param[distribution][0]
            self.loc = fit.fitted_param[distribution][1]
            self.scale = fit.fitted_param[distribution][2]
        elif distribution == "exponpow":
            self.distribution = "exponpow"
            self.b = fit.fitted_param[distribution][0]
            self.loc = fit.fitted_param[distribution][1]
            self.scale = fit.fitted_param[distribution][2]
        elif distribution == "rayleigh":
            self.distribution = "rayleigh"
            self.loc = fit.fitted_param[distribution][0]
            self.scale = fit.fitted_param[distribution][1]
        return self
    
    def precision(self):
        """
        A technique to provide the precision of the ThunderSTORM values, determined by the best fit of the histogram.

        Return:
            A plot displaying the histogram and important variables.
        """
        fig0, ax0 = plt.subplots(figsize=(10,5))
        ax0.set_title("Precision Analysis")
        ax0.set_ylabel("Density")
        ax0.set_xlabel(self.data)
        sns.histplot(x=self.points[self.data], 
                     ax=ax0, stat="density", 
                     kde=True, bins = 100, 
                     linewidth=0.5, label=("N="+str(len(self.points[self.data]))), line_kws={'label': "Kernel Density Estimation (KDE) Fit"})
        if self.distribution == "norm":
            fitting = norm.pdf(self.Xaxis, loc=self.loc, scale=self.scale)
        elif self.distribution == "lognorm":
            fitting = lognorm.pdf(self.Xaxis, s=self.s, loc=self.loc, scale=self.scale)
        elif self.distribution == "gamma":
            fitting = gamma.pdf(self.Xaxis, a=self.a, loc=self.loc, scale=self.scale)
        elif self.distribution == "chi2":
            fitting = chi2.pdf(self.Xaxis, df=self.df, loc=self.loc, scale=self.scale)
        elif self.distribution == "exponpow":
            fitting = exponpow.pdf(self.Xaxis, b=self.b, loc=self.loc, scale=self.scale)
        elif self.distribution == "rayleigh":
            fitting = rayleigh.pdf(self.Xaxis, loc=self.loc, scale=self.scale)
        ax0.plot(self.Xaxis, fitting, color="red", linestyle="dashed", label="Probability Density Function (PDF) Fit")
        a = round(self.Xaxis[np.where(fitting == fitting.max())][0], 4)
        b = round(fitting[np.where(fitting == fitting.max())][0], 4)
        ax0.scatter(a, b, color="red", label=("["+str(a)+", "+str(b)+"]"))
        ax0.legend()
        plt.show()
        return self
