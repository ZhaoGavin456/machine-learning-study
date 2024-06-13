import numpy as np
import matplotlib.pyplot as plt
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, confusion_matrix, classification_report

# Step 1: Generate synthetic data
X, y = make_classification(n_samples=1000, n_features=2, n_informative=2, n_redundant=0,
                           n_clusters_per_class=1, random_state=42)

# Add intercept term to X
X = np.hstack((np.ones((X.shape[0], 1)), X))

# Step 2: Split the data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

# Step 3: Define helper functions
def sigmoid(z):
    return 1 / (1 + np.exp(-z))

def cost_function(X, y, weights):
    m = X.shape[0]
    h = sigmoid(np.dot(X, weights))
    return (-1/m) * np.sum(y * np.log(h) + (1 - y) * np.log(1 - h))

def gradient_ascent(X, y, weights, learning_rate, iterations):
    m = X.shape[0]
    cost_history = []

    for i in range(iterations):
        h = sigmoid(np.dot(X, weights))
        gradient = np.dot(X.T, (y - h)) / m
        weights += learning_rate * gradient
        cost = cost_function(X, y, weights)
        cost_history.append(cost)

    return weights, cost_history

# Step 4: Train the model using batch gradient ascent
learning_rate = 0.1
iterations = 1000
weights = np.zeros(X_train.shape[1])
weights, cost_history = gradient_ascent(X_train, y_train, weights, learning_rate, iterations)

# Step 5: Make predictions
def predict(X, weights):
    return np.round(sigmoid(np.dot(X, weights)))

y_pred = predict(X_test, weights)

# Step 6: Evaluate the model
accuracy = accuracy_score(y_test, y_pred)
conf_matrix = confusion_matrix(y_test, y_pred)
class_report = classification_report(y_test, y_pred)

print(f"Accuracy: {accuracy}")
print("Confusion Matrix:")
print(conf_matrix)
print("Classification Report:")
print(class_report)

# Step 7: Plot the cost function history and decision boundary
fig, ax = plt.subplots(1, 2, figsize=(14, 6))

# Plot the cost function history
ax[0].plot(range(iterations), cost_history)
ax[0].set_xlabel("Iterations")
ax[0].set_ylabel("Cost")
ax[0].set_title("Cost Function History")
ax[0].grid(True)

# Plot the decision boundary with training and test points
def plot_decision_boundary(ax, X_train, y_train, X_test, y_test, weights):
    x_min, x_max = X_train[:, 1].min() - 1, X_train[:, 1].max() + 1
    y_min, y_max = X_train[:, 2].min() - 1, X_train[:, 2].max() + 1
    xx, yy = np.meshgrid(np.arange(x_min, x_max, 0.01),
                         np.arange(y_min, y_max, 0.01))
    Z = predict(np.c_[np.ones((xx.ravel().shape[0], 1)), xx.ravel(), yy.ravel()], weights)
    Z = Z.reshape(xx.shape)
    ax.contourf(xx, yy, Z, alpha=0.3, cmap=plt.cm.Paired)

    # Plot training points
    scatter1 = ax.scatter(X_train[:, 1], X_train[:, 2], c=y_train, marker='o', edgecolor='k', s=50, cmap=plt.cm.Paired, label='Train')
    # Plot test points
    scatter2 = ax.scatter(X_test[:, 1], X_test[:, 2], c=y_test, marker='x', s=80, cmap=plt.cm.Paired, label='Test')

    ax.set_xlim(xx.min(), xx.max())
    ax.set_ylim(yy.min(), yy.max())
    ax.set_title("Logistic Regression Decision Boundary")
    ax.set_xlabel("Feature 1")
    ax.set_ylabel("Feature 2")
    ax.legend()
    ax.grid(True)

plot_decision_boundary(ax[1], X_train, y_train, X_test, y_test, weights)

plt.tight_layout()
plt.show()
