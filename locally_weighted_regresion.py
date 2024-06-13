import numpy as np
import matplotlib.pyplot as plt


def local_weighted_regression(x, y, x0, tau):
    # Calculate weights
    weights = np.exp(-((x - x0) ** 2) / (2 * tau ** 2))

    # Expand dimensions for matrix operations
    x = x[:, np.newaxis]
    y = y[:, np.newaxis]
    weights = weights[:, np.newaxis]

    # Design matrix
    X = np.hstack((np.ones((len(x), 1)), x))

    # Compute beta
    beta = np.linalg.pinv(X.T @ (weights * X)) @ (X.T @ (weights * y))

    return beta[0] + beta[1] * x0


def lowess(x, y, tau, x_vals=None):
    if x_vals is None:
        x_vals = np.linspace(min(x), max(x), len(x))
    y_vals = np.array([local_weighted_regression(x, y, x0, tau) for x0 in x_vals])
    return x_vals, y_vals


# Example data
x = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
y = np.array([2.5, 2.0, 3.5, 4.0, 6.0, 6.0, 7.5, 7.0, 8.0, 9.0])

# Bandwidth parameter
tau = 1.0

# Perform LOWESS
x_vals, y_vals = lowess(x, y, tau)

# Test points
test_points = [3.5, 7.5]
test_point_y = [local_weighted_regression(x, y, test_point, tau) for test_point in test_points]
# test_point_y_str = [f'{y:.2f}' for y in test_point_y]  # Convert to string
test_point_y_str = [str(y) for y in test_point_y]  # Convert to string

# Plot the results
plt.scatter(x, y, color='red', label='Data')
plt.plot(x_vals, y_vals, color='blue', label='LOWESS')
plt.scatter(test_points, test_point_y, color='green', label='Test Points')
for i, txt in enumerate(test_point_y_str):
    plt.annotate(txt, (test_points[i], test_point_y[i]), textcoords="offset points", xytext=(0,10), ha='center')
plt.xlabel('X')
plt.ylabel('Y')
plt.legend()
plt.show()
