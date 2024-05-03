

import numpy as np
from sklearn.covariance import MinCovDet
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

""" def generate_3d_data(n_samples):
    """
    Generate random 3D data
    """
    np.random.seed(0)
    data = np.random.randn(n_samples, 3)
    return data

def fit_3d_mce(data):
    """
    Fit Minimum Covariance Ellipsoid to 3D data
    """
    mcd = MinCovDet()
    mcd.fit(data)
    return mcd.location_, mcd.covariance_

def plot_3d_mce(data, center, cov):
    """
    Plot 3D data and Minimum Covariance Ellipsoid
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], c='b', marker='o', label='Data points')
    
    # Plot ellipsoid
    u, s, vh = np.linalg.svd(cov)
    radii = np.sqrt(s)
    u *= radii[0]
    vh *= radii[1]
    vh = vh.T
    for i in range(0, 360, 45):
        phi = i * np.pi / 180
        r = np.array([[np.cos(phi), -np.sin(phi)], [np.sin(phi), np.cos(phi)]])
        ellipse = np.dot(u, np.dot(r, vh)) + center
        ax.plot(ellipse[:, 0], ellipse[:, 1], ellipse[:, 2], 'r')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.legend()
    plt.show()

# Generate 3D data
data = generate_3d_data(100)

print(data)
# Fit MCE to the data
center, cov = fit_3d_mce(data)
print(center)
# Plot data and MCE
plot_3d_mce(data, center, cov) """








