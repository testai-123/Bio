import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.metrics import accuracy_score, classification_report
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
data = pd.read_csv('Genomics.csv')

# Select features and target
features = data[['area', 'area_type', 'specimens', 'percentage', 'specimens_7d_avg', 'percentage_7d_avg']]
labels = data['variant_name']

# Encode categorical features
encoder = LabelEncoder()
features['area'] = encoder.fit_transform(features['area'])
features['area_type'] = encoder.fit_transform(features['area_type'])

# Handle missing values
imputer = SimpleImputer(strategy='mean')
features = imputer.fit_transform(features)

# Split data
X_train, X_test, y_train, y_test = train_test_split(features, labels, test_size=0.2, random_state=42)

# Scale features
scaler = StandardScaler()
X_train_scaled = scaler.fit_transform(X_train)
X_test_scaled = scaler.transform(X_test)

# Train Random Forest Classifier
rf_classifier = RandomForestClassifier(n_estimators=100, max_depth=None, random_state=42)
rf_classifier.fit(X_train_scaled, y_train)

# Evaluate Random Forest
rf_predictions = rf_classifier.predict(X_test_scaled)
print("Random Forest Accuracy:", accuracy_score(y_test, rf_predictions))
print("Random Forest Classification Report:\n", classification_report(y_test, rf_predictions))

# Plot top 5 feature importances
rf_importances = rf_classifier.feature_importances_
feature_names = ['area', 'area_type', 'specimens', 'percentage', 'specimens_7d_avg', 'percentage_7d_avg']
feature_importance_df = pd.DataFrame({'Feature': feature_names, 'Importance': rf_importances}).sort_values(by='Importance', ascending=False)
sns.barplot(x='Importance', y='Feature', data=feature_importance_df.head(5))
plt.title('Top 5 Important Features (Random Forest)')
plt.show()

# Train Support Vector Machine Classifier
svm_model = SVC(C=1.0, kernel='rbf', gamma='scale', probability=True, random_state=42)
svm_model.fit(X_train_scaled, y_train)

# Evaluate SVM
svm_predictions = svm_model.predict(X_test_scaled)
print("SVM Accuracy:", accuracy_score(y_test, svm_predictions))
print("SVM Classification Report:\n", classification_report(y_test, svm_predictions))


# Fin
print("Random Forest Accuracy:", accuracy_score(y_test, rf_predictions))
print("SVM Accuracy:", accuracy_score(y_test, svm_predictions))

"""1. Random Forest Classifier:
  * Accuracy: 41.98%
  * Precision: Higher for some classes (e.g., "Alpha" at 0.62, "Omicron" at 0.72) but inconsistent across others.
  * Recall: Highest for "Beta" (0.83) and "Other" (0.81), indicating it detected these variants better.
  * F1-Score: Low for most classes, indicating imbalanced prediction quality.
  * Strengths: Performed reasonably on "Other" and "Omicron" classes.
  * Weaknesses: Struggled with classes like "Lambda" and "Mu".
2. Support Vector Machine (SVM):
  * Accuracy: 23.49%
  * Precision: Very low for most classes, with "Alpha" being the exception (0.83).
  * Recall: Significantly high only for "Lambda" (0.99) and "Total" (0.98), showing a bias towards these classes.
  * F1-Score: Very low across most categories, indicating poor prediction quality.
  * Strengths: Only high recall for "Lambda" and "Total", indicating some true positives were detected.
  * Weaknesses: Very poor performance in predicting most other classes.
"""

