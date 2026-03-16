from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.svm import SVC

model_params = {
    "lr": {
        "estimator": LogisticRegression(max_iter=300, class_weight="balanced"),
    },
    "svm": {
        "estimator": SVC(kernel="linear", class_weight="balanced"),
        "param_grid": {
            "ensemble": {
                "base_estimator__estimator__C": [0.001, 0.05, 0.01, 0.05, 0.1, 0.5, 1],
                "subset_frac": [0.5, 0.75, 1.0]
            },
            "nonensemble": {
                "estimator__C": [0.001, 0.05, 0.01, 0.05, 0.1, 0.5, 1],
            }
        }
    },
    "rf": {
        "estimator": RandomForestClassifier(class_weight="balanced"),
        "param_grid": {
            # we label rf as "nonensemble" since it is not using our random undersampling 
            "nonensemble": {"estimator__n_estimators": [100, 300, 500], "estimator__max_depth": [3, 4, 5, 6, 7]}
        }
    },
}