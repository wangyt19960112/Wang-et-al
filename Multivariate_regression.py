"""
Multivariate linear regression analysis of signaling pathway activity
and lineage composition in axial organoid models.

Input:
    data/Signal_protocol.xlsx

Output:
    results/regression/
        *_fit_colored.svg
        lineage_signal_weights_with_intercept.svg
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score, mean_absolute_error
from sklearn.preprocessing import StandardScaler


# ============================================================
# 1. File paths
# ============================================================

DATA_DIR = Path("data")
RESULTS_DIR = Path("results") / "regression"

INPUT_FILE = DATA_DIR / "Signal_protocol.xlsx"

# Create output directory if it does not exist
RESULTS_DIR.mkdir(parents=True, exist_ok=True)


# ============================================================
# 2. Load data
# ============================================================

df = pd.read_excel(INPUT_FILE)

print("Columns in input file:")
print(df.columns.tolist())


# Combine hindbrain and cervical-thoracic neural tube
df["Posterior NT(%)"] = (
    df["Hindbrain(%)"].fillna(0)
    + df["Cervical-Thoracic NT(%)"].fillna(0)
)


# ============================================================
# 3. Calculate signaling pathway activity
# ============================================================

def zscore(series):
    """
    Calculate z-score using population standard deviation (ddof=0).
    """
    return (series - series.mean()) / series.std(ddof=0)


# WNT signaling
df["WNT_signal"] = zscore(
    df["CHIR (μM)"] * df["CHIR_time"]
)


# FGF signaling
df["FGF2_z"] = zscore(
    df["FGF2 (ng ml-1)"] * df["FGF2_time"]
)

df["FGF8_z"] = zscore(
    df["FGF8 (ng ml-1)"] * df["FGF8_time"]
)

df["FGF_signal"] = (
    df["FGF2_z"]
    + df["FGF8_z"]
)


# Retinoic acid signaling
df["RA_z"] = zscore(
    df["RA (nM)"] * df["RA_time"]
)

df["RAL_z"] = zscore(
    df["RAL(μM)"] * df["RAL_time"]
)

df["RA_signal"] = (
    df["RA_z"]
    + df["RAL_z"]
)


# Ventral / SHH signaling
# Calculated here for completeness but not included
# in the regression model below.
df["Neural_ventral"] = zscore(
    df["SAG(μM)"] * df["SAG_time"]
)


# SMAD inhibition
smadi_cols = [
    "SMADi_LDN (μM)",
    "SMADi_SB (μM)",
    "SMADi_DMH1 (μM)",
    "SMADi_A8301(μM)",
]

smadi_time_cols = [
    "SMADi_LDN_time",
    "SMADi_SB_time",
    "SMADi_DMH1_time",
    "SMADi_A8301_time",
]

df["SMADi_signal"] = 0.0

for conc_col, time_col in zip(smadi_cols, smadi_time_cols):

    if conc_col in df.columns and time_col in df.columns:

        exposure = df[conc_col] * df[time_col]

        df["SMADi_signal"] += zscore(exposure)


# ============================================================
# 4. Define regression features and lineage outputs
# ============================================================

features = [
    "WNT_signal",
    "FGF_signal",
    "RA_signal",
    "SMADi_signal",
]

lineage_cols = [
    "Anterior NT(%)",
    "Posterior NT(%)",
    "NMP(%)",
    "Mesoderm(%)",
]


X = df[features]


# Standardize pathway-level features before regression
scaler = StandardScaler()

X_scaled = scaler.fit_transform(X)


# ============================================================
# 5. Multivariate linear regression
# ============================================================

models = {}
coef_dict = {}
predictions = {}


for lineage in lineage_cols:

    y = df[lineage]

    model = LinearRegression(
        fit_intercept=True
    )

    model.fit(
        X_scaled,
        y
    )

    y_pred = model.predict(
        X_scaled
    )


    # Store model and predictions
    models[lineage] = model
    predictions[lineage] = y_pred


    lineage_name = lineage.replace("(%)", "").strip()

    coef_dict[lineage_name] = model.coef_

    df[f"{lineage_name}_fit"] = y_pred


    # --------------------------------------------------------
    # Print regression equation and performance
    # --------------------------------------------------------

    equation_terms = [
        f"{model.intercept_:.3f}"
    ]

    for feature, coefficient in zip(
        features,
        model.coef_
    ):

        equation_terms.append(
            f"{coefficient:.3f} × {feature}"
        )


    print(f"\n=== {lineage} ===")

    print(
        "Lineage = "
        + " + ".join(equation_terms)
    )

    print(
        f"R² = {r2_score(y, y_pred):.2f}, "
        f"MAE = {mean_absolute_error(y, y_pred):.3f}"
    )


# ============================================================
# 6. Plot predicted versus observed lineage proportions
# ============================================================

# Keep text editable in SVG files
plt.rcParams["svg.fonttype"] = "none"


def get_color(sample_id):
    """
    Assign plotting colors according to Sample ID prefix.
    """

    sample_id = str(sample_id)

    if sample_id.startswith("SM"):
        return "green"

    elif sample_id.startswith("NT"):
        return "red"

    elif sample_id.startswith("TK"):
        return "orange"

    else:
        return "grey"


colors = df["Sample ID"].apply(get_color)


for lineage in lineage_cols:

    lineage_name = lineage.replace("(%)", "").strip()

    y_true = df[lineage]

    y_pred = df[f"{lineage_name}_fit"]


    fig, ax = plt.subplots(
        figsize=(4, 4)
    )


    # Observed versus predicted values
    ax.scatter(
        y_true,
        y_pred,
        s=70,
        c=colors
    )


    # Identity line (y = x)
    min_v = min(
        y_true.min(),
        y_pred.min()
    )

    max_v = max(
        y_true.max(),
        y_pred.max()
    )

    ax.plot(
        [min_v, max_v],
        [min_v, max_v],
        "k--"
    )


    # Sample labels
    for i, sample_id in enumerate(df["Sample ID"]):

        ax.text(
            y_true.iloc[i] + 1,
            y_pred.iloc[i],
            str(sample_id),
            fontsize=8
        )


    ax.set_xlabel(
        "Observed (%)"
    )

    ax.set_ylabel(
        "Predicted (%)"
    )

    ax.set_title(
        lineage_name
    )


    fig.tight_layout()


    output_name = (
        lineage_name
        .replace(" ", "_")
        + "_fit_colored.svg"
    )


    fig.savefig(
        RESULTS_DIR / output_name,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close(fig)


# ============================================================
# 7. Plot regression coefficients
# ============================================================

coef_df = pd.DataFrame(
    coef_dict,
    index=features
)


fig, ax = plt.subplots(
    figsize=(8, 5)
)


sns.heatmap(
    coef_df.T,
    cmap="coolwarm",
    center=0,
    annot=True,
    fmt=".2f",
    ax=ax
)


ax.set_xlabel(
    "Signaling pathway"
)

ax.set_ylabel(
    "Lineage"
)

ax.set_title(
    "Linear regression coefficients"
)


fig.tight_layout()


fig.savefig(
    RESULTS_DIR / "lineage_signal_weights_with_intercept.svg",
    dpi=300,
    bbox_inches="tight"
)

plt.close(fig)


# ============================================================
# 8. Save regression results
# ============================================================

df.to_csv(
    RESULTS_DIR / "regression_results.csv",
    index=False
)

coef_df.T.to_csv(
    RESULTS_DIR / "regression_coefficients.csv"
)


print(
    f"\nAnalysis completed. Results saved to: {RESULTS_DIR}"
)