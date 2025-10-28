import gradio as gr
import pandas as pd
import joblib
import numpy as np
from sksurv.nonparametric import kaplan_meier_estimator
import warnings

# --- Global Settings & Model Loading ---
warnings.filterwarnings("ignore")

# Define the median risk score from your validation script
# This is the threshold to separate High/Low risk groups
#
MEDIAN_RISK_THRESHOLD = -0.0452

def load_model_and_features():
    """
    Loads the trained model and extracts key information.
   
    """
    try:
        # Load the saved LASSO-Cox model
        #
        model = joblib.load("trained_lasso_cox_model.pkl")
        
        # Get the 720 features the model was trained on
        # [cite: 21-68]
        required_features = model.feature_names_in_
        
        # Get the 44-gene signature (features with non-zero coefficients)
        #
        coef = pd.Series(model.coef_, index=required_features)
        signature_genes = coef[coef != 0].sort_values(ascending=False)
        
        signature_df = pd.DataFrame({
            "Gene_ID": signature_genes.index,
            "LASSO_Coefficient": signature_genes.values
        })
        
        return model, required_features, signature_df
    except FileNotFoundError:
        print("ERROR: Model file 'trained_lasso_cox_model.pkl' not found.")
        return None, [], pd.DataFrame()
    except Exception as e:
        print(f"Error loading model: {e}")
        return None, [], pd.DataFrame()

# Load model and data on app startup
MODEL, REQUIRED_FEATURES, SIGNATURE_DF = load_model_and_features()
EXAMPLE_FILE = "example_patient_data.csv"

# --- Core Prediction Function ---

def predict_risk(patient_file):
    """
    Takes an uploaded CSV, validates it, and returns risk scores.
    """
    if MODEL is None or not hasattr(patient_file, 'name'):
        raise gr.Error("Model is not loaded. Please check server logs.")
        
    try:
        # Read the uploaded CSV file
        df = pd.read_csv(patient_file.name)
    except Exception as e:
        raise gr.Error(f"Failed to read CSV: {e}")

    # --- 1. Validation ---
    # Check for a Patient ID column
    id_col = next((col for col in df.columns if 'patient' in col.lower() or 'id' in col.lower()), None)
    if id_col is None:
        gr.Warning("No 'Patient_ID' column found. Using row index as ID.")
        df.insert(0, "Patient_ID", [f"Patient_{i+1}" for i in df.index])
        id_col = "Patient_ID"
        
    # Check if all 720 required genes are in the file
    missing_cols = set(REQUIRED_FEATURES) - set(df.columns)
    if missing_cols:
        raise gr.Error(
            f"Input file is missing {len(missing_cols)} required gene(s). "
            f"First missing gene: {list(missing_cols)[0]}. "
            "Please download the example file for the correct format."
        )

    # --- 2. Pre-processing & Prediction ---
    try:
        # Ensure data is in the correct order (all 720 features)
        X_new = df[REQUIRED_FEATURES]
        
        # Scale the data using the *exact same* scaler from training
        #
        X_scaled = MODEL.scaler.transform(X_new)
        
        # Predict the prognostic risk score
        #
        risk_scores = MODEL.predict(X_scaled)
        
    except Exception as e:
        raise gr.Error(f"Error during model prediction: {e}")

    # --- 3. Format Output ---
    results_df = pd.DataFrame({
        "Patient ID": df[id_col],
        "Prognostic Risk Score": np.round(risk_scores, 4),
        "Risk Group": ["High Risk" if score > MEDIAN_RISK_THRESHOLD else "Low Risk" for score in risk_scores]
    })
    
    # Sort by risk score (highest risk first)
    results_df = results_df.sort_values(by="Prognostic Risk Score", ascending=False)
    
    return results_df

# --- Gradio UI ---

# Using Blocks for a more professional, multi-tab layout
with gr.Blocks(theme=gr.themes.Soft(), title="CRC Prognostic Calculator") as app:
    
    gr.Markdown(
        """
        # 🧬 CRC-RBM47 Prognostic Signature Calculator
        This tool predicts patient prognosis in Colorectal Cancer (CRC) based on a machine learning model.
        
        The model is a **LASSO-Cox Proportional Hazards model** [cite: 14] trained on TCGA-COAD patient data[cite: 13]. 
        It uses a 44-gene signature derived from the **RBM47 regulatory network**, first identified in an *in-vitro* study[cite: 5, 7, 11].
        
        **Project by:** Abhishek S R, MSc Bioinformatics [cite: 4]
        **Source Project:** "To develop an Integrated Pipeline to Validate and Model the RBM47-Regulated Gene signature and Splicing Signature in Colorectal Cancer" 
        """
    )
    
    with gr.Tabs():
        
        # --- TAB 1: Run Prediction ---
        with gr.TabItem("Run Prognostic Prediction", id=0):
            gr.Markdown("Upload a CSV file containing normalized gene expression data for one or more patients. The file must contain columns for all 720 genes used in training.")
            
            with gr.Row():
                with gr.Column(scale=1):
                    file_input = gr.File(label="Upload Patient Expression CSV", file_types=[".csv"])
                    
                    gr.File(
                        label="Download Example File (Required Format)",
                        value=EXAMPLE_FILE,
                        interactive=False
                    )
                    
                    predict_btn = gr.Button("Calculate Risk Score", variant="primary")
                    
                with gr.Column(scale=2):
                    gr.Markdown("### Prediction Results")
                    output_df = gr.DataFrame(
                        label="Patient Risk Scores",
                        headers=["Patient ID", "Prognostic Risk Score", "Risk Group"],
                        datatype=["str", "number", "str"]
                    )
            
            # Wire up the prediction
            predict_btn.click(
                fn=predict_risk,
                inputs=[file_input],
                outputs=[output_df]
            )

        # --- TAB 2: Model Details & Validation ---
        with gr.TabItem("About the Model & Validation", id=1):
            gr.Markdown("### Model Performance & Validation")
            gr.Markdown(
                f"""
                This model was trained and validated on **410** TCGA-COAD patients.
                
                * **Model Type:** LASSO-Cox Proportional Hazards [cite: 14]
                * **Model Performance (Test C-index):** **0.7314**
                * **Validation Log-Rank P-value:** **9.71e-17** (p < 0.001)
                """
            )
            
            # Display the "Money Plot"
            gr.Image(
                "Validation plot of ML model.png", 
                label="Prognostic Signature Validation (n=410)", 
                interactive=False
            )
            
            gr.Markdown(
                """
                ### 44-Gene Prognostic Signature
                The final LASSO-Cox model selected the following 44 genes (out of 720 candidates) as having predictive prognostic value.
                A positive coefficient contributes to a **'High Risk'** score, while a negative coefficient contributes to a **'Low Risk'** score.
               
                """
            )
            
            gr.DataFrame(
                value=SIGNATURE_DF,
                label="Final 44-Gene Signature",
                interactive=False
            )

# Launch the app
app.launch()