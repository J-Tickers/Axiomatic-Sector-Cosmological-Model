import numpy as np
from cobaya.likelihood import Likelihood
from scipy.integrate import trapezoid
from scipy.special import eval_legendre
from scipy.interpolate import interp1d

class AuroraPk(Likelihood):
    def initialize(self):
        DATA_PATH = "new_datasets/des_y6/y6_bao/Data/Clustering/DESY6BAO_datavectors/acf/acf_data_bin1.dat"
        COV_PATH = "new_datasets/des_y6/y6_bao/Covariances/acf_cov_Planck_datalike.txt"
        NZ_PATH = "new_datasets/des_y6/y6_bao/Data/Redshift_distributions/nz_fiducial.txt"

        self.data_array = np.loadtxt(DATA_PATH)
        self.theta_deg = self.data_array[:, 0]
        self.data_vector = self.data_array[:, 1]
        cov_matrix = np.loadtxt(COV_PATH)

        n_bins = len(self.data_vector)
        self.inv_cov = np.linalg.inv(cov_matrix[:n_bins, :n_bins])

        nz_data = np.loadtxt(NZ_PATH)
        self.nz_interp = interp1d(nz_data[:, 0], nz_data[:, 1], bounds_error=False, fill_value=0.0)

    def get_requirements(self):
        return {
            "angular_diameter_distance": {"z": np.linspace(0.1, 1.5, 50)},
            "Hubble": {"z": np.linspace(0.1, 1.5, 50)},
            "Pk_interpolator": {"z": np.linspace(0.1, 1.5, 50), "k_max": 5.0, "nonlinear": True}
        }

    def logp(self, **params_values):
        provider = self.provider

        z_arr = np.linspace(0.1, 1.5, 50)
        chi_arr = np.array([provider.get_angular_diameter_distance(z) * (1+z) for z in z_arr])
        hz_arr = np.array([provider.get_Hubble(z) for z in z_arr])
        nz_arr = self.nz_interp(z_arr)

        pk_interp = provider.get_Pk_interpolator(nonlinear=True)

        ell_arr = np.arange(10, 1000, 20)
        c_ell = np.zeros_like(ell_arr, dtype=float)

        for i, ell in enumerate(ell_arr):
            integrand = np.zeros_like(z_arr)
            for j, z in enumerate(z_arr):
                k = (ell + 0.5) / chi_arr[j]
                pk = pk_interp.P(z, k)
                integrand[j] = (nz_arr[j] / chi_arr[j])**2 * pk / hz_arr[j]
            c_ell[i] = trapezoid(integrand, x=z_arr)

        theta_rad = np.radians(self.theta_deg)
        w_theta = np.zeros(len(theta_rad))

        for i, t in enumerate(theta_rad):
            p_ell = eval_legendre(ell_arr, np.cos(t))
            w_theta[i] = np.sum((2 * ell_arr + 1) / (4 * np.pi) * p_ell * c_ell)

        residual = self.data_vector - w_theta
        chi2 = np.dot(residual, np.dot(self.inv_cov, residual))
        
        return -0.5 * chi2
