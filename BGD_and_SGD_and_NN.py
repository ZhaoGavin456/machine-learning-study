import numpy as np
import matplotlib.pyplot as plt
import tensorflow as tf
from sklearn.preprocessing import StandardScaler

def batch_gradient_descent(X, y, learning_rate=0.01, iterations=1000):
    m, n = X.shape
    X = np.c_[np.ones((m, 1)), X]  # Add bias term
    theta = np.zeros(n + 1)
    for iteration in range(iterations):
        gradients = 2/m * X.T.dot(X.dot(theta) - y)
        theta -= learning_rate * gradients
    return theta

def stochastic_gradient_descent(X, y, learning_rate=0.01, iterations=1000):
    m, n = X.shape
    X = np.c_[np.ones((m, 1)), X]  # Add bias term
    theta = np.zeros(n + 1)
    for iteration in range(iterations):
        for i in range(m):
            random_index = np.random.randint(m)
            xi = X[random_index:random_index+1]
            yi = y[random_index:random_index+1]
            gradients = 2 * xi.T.dot(xi.dot(theta) - yi)
            theta -= learning_rate * gradients
    return theta

def train_linear_nn(X, y, learning_rate=0.01, epochs=500):
    model = tf.keras.Sequential([
        tf.keras.layers.Dense(1, input_shape=(1,), activation='linear')
    ])
    model.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate), loss='mse')
    model.fit(X, y, epochs=epochs, verbose=0)
    return model

def plot_regression_lines(X, y, y_pred_bgd, y_pred_sgd, y_pred_nn, y_pred_linear_nn):
    plt.scatter(X, y, color='blue', label='Training data')
    plt.plot(X, y_pred_bgd, color='red', label='BGD Regression line')
    plt.plot(X, y_pred_sgd, color='green', label='SGD Regression line')
    plt.plot(X, y_pred_nn, color='purple', label='NN Regression line')
    plt.plot(X, y_pred_linear_nn, color='orange', label='Linear NN Regression line')
    plt.xlabel('X')
    plt.ylabel('y')
    plt.legend()
    plt.show()

def mean_squared_error(y_true, y_pred):
    mse = np.mean((y_true - y_pred) ** 2)
    return mse

# Generate synthetic data
np.random.seed(42)  # For reproducibility
X = np.linspace(1, 50, 50).reshape(-1, 1)
y = 3 * X.squeeze() + np.random.randn(50) * 10  # Linear relationship with noise

# Normalize the data
scaler_X = StandardScaler()
scaler_y = StandardScaler()
X_normalized = scaler_X.fit_transform(X)
y_normalized = scaler_y.fit_transform(y.reshape(-1, 1)).flatten()

# BGD
theta_bgd = batch_gradient_descent(X_normalized, y_normalized, learning_rate=0.01, iterations=10000)
y_pred_bgd_normalized = np.c_[np.ones((X_normalized.shape[0], 1)), X_normalized].dot(theta_bgd)
y_pred_bgd = scaler_y.inverse_transform(y_pred_bgd_normalized.reshape(-1, 1)).flatten()
mse_bgd = mean_squared_error(y, y_pred_bgd)
print("BGD Theta:", theta_bgd)
print("BGD Mean Squared Error:", mse_bgd)

# SGD
theta_sgd = stochastic_gradient_descent(X_normalized, y_normalized, learning_rate=0.01, iterations=10000)
y_pred_sgd_normalized = np.c_[np.ones((X_normalized.shape[0], 1)), X_normalized].dot(theta_sgd)
y_pred_sgd = scaler_y.inverse_transform(y_pred_sgd_normalized.reshape(-1, 1)).flatten()
mse_sgd = mean_squared_error(y, y_pred_sgd)
print("SGD Theta:", theta_sgd)
print("SGD Mean Squared Error:", mse_sgd)

# Build and train the neural network model (with hidden layers)
model_nn = tf.keras.Sequential([
    tf.keras.layers.Dense(64, input_shape=(1,), activation='relu'),
    tf.keras.layers.Dense(32, activation='relu'),
    tf.keras.layers.Dense(1, activation='linear')
])
model_nn.compile(optimizer=tf.keras.optimizers.Adam(learning_rate=0.01), loss='mse')
model_nn.fit(X_normalized, y_normalized, epochs=500, verbose=0)
y_pred_nn_normalized = model_nn.predict(X_normalized)
y_pred_nn = scaler_y.inverse_transform(y_pred_nn_normalized).flatten()
mse_nn = model_nn.evaluate(X_normalized, y_normalized, verbose=0)
print(f'NN Mean Squared Error: {mse_nn}')

# Build and train the linear neural network model (without hidden layers)
model_linear_nn = train_linear_nn(X_normalized, y_normalized, learning_rate=0.01, epochs=500)
y_pred_linear_nn_normalized = model_linear_nn.predict(X_normalized)
y_pred_linear_nn = scaler_y.inverse_transform(y_pred_linear_nn_normalized).flatten()
mse_linear_nn = model_linear_nn.evaluate(X_normalized, y_normalized, verbose=0)
print(f'Linear NN Mean Squared Error: {mse_linear_nn}')

# Plot the results
plot_regression_lines(X, y, y_pred_bgd, y_pred_sgd, y_pred_nn, y_pred_linear_nn)
