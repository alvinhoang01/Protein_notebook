import streamlit as st
import pandas as pd
import numpy as np
import os
import plotly.express as px

def show_protein_total_viewer():
    # --- Data Loading and Processing ---
    
    # Load the new data file
    data_path = os.path.join('data', 'Protein_total_for website.csv')
    data = pd.read_csv(data_path)

    # Pivot the data into long format
    data_long = pd.melt(
        data, 
        id_vars=['PG.ProteinAccessions'], 
        value_vars=[
            'Plasma_R1', 'Plasma_R2', 'Plasma_R3',
            'Serum_R1', 'Serum_R2', 'Serum_R3',
            'Whole.blood_R1', 'Whole.blood_R2', 'Whole.blood_R3',
            'Tissue_R1', 'Tissue_R2', 'Tissue_R3'
        ], 
        var_name='Sample', 
        value_name='Intensity'
    )

    # Create 'Matrice' and 'Replicate' columns based on patterns in 'Sample'
    data_long['Matrice'] = np.select(
        [
            data_long['Sample'].str.contains('Plasma'),
            data_long['Sample'].str.contains('Serum'),
            data_long['Sample'].str.contains('Whole.blood'),
            data_long['Sample'].str.contains('Tissue')
        ],
        ['Plasma', 'Serum', 'Whole Blood', 'Tissue'],
        default=np.nan
    )

    data_long['Replicate'] = np.select(
        [
            data_long['Sample'].str.contains('R1'),
            data_long['Sample'].str.contains('R2'),
            data_long['Sample'].str.contains('R3')
        ],
        ['R1', 'R2', 'R3'],
        default=np.nan
    )

    # Calculate mean, standard deviation, and CV for each Protein-Matrice group
    data_cv = data_long.groupby(['PG.ProteinAccessions', 'Matrice']).agg(
        mean_value=('Intensity', 'mean'),
        sd_value=('Intensity', 'std')
    ).reset_index()
    data_cv['CV'] = (data_cv['sd_value'] / data_cv['mean_value']) * 100

    # Merge the CV back into the long data frame
    data_long = pd.merge(
        data_long, 
        data_cv[['PG.ProteinAccessions', 'Matrice', 'CV']], 
        on=['PG.ProteinAccessions', 'Matrice'], 
        how='left'
    )

    # --- Streamlit App Interface ---
    st.title("Proteomic Matrices")

    # Sidebar inputs for protein search and selection
    st.sidebar.header("Search Options")
    protein_options = data_long['PG.ProteinAccessions'].unique()
    protein_search = st.sidebar.text_input("Enter PG.ProteinAccessions", key="filtered_protein_search")
    protein_select = st.sidebar.selectbox("Select PG.ProteinAccessions", protein_options, key="filtered_protein_select")
    search = st.sidebar.button("Search", key="filtered_search_button")
    reset = st.sidebar.button("Reset to Full Table", key="filtered_reset_button")

    # Display either full or filtered data based on search or reset
    if search:
        # Determine which protein to filter by
        protein = protein_search if protein_search != "" else protein_select
        filtered_data = data_long[data_long['PG.ProteinAccessions'] == protein]
        filtered_data_cv = data_cv[data_cv['PG.ProteinAccessions'] == protein]

        # Display filtered data and plots
        st.subheader(f"Data Table for {protein}")
        st.dataframe(filtered_data, use_container_width=True)  # Use full width for filtered data

        # Create tabs for the two plots
        tab1, tab2 = st.tabs(["MS Intensity Plot", "CV Plot"])

        with tab1:
            # Create the intensity plot with adjusted settings
            intensity_fig = px.bar(
                filtered_data,
                x='Matrice',
                y='Intensity',
                color='Matrice',
                hover_data=['Replicate'],
                labels={'Intensity': 'MS Intensity'},
                title=f"MS Intensity for {protein}"
            )
            
            # Update layout for alignment - removes grouped spacing
            intensity_fig.update_layout(
                barmode='overlay',  # Center bars over each label
                xaxis=dict(
                    categoryorder="array",
                    categoryarray=["Plasma", "Serum", "Whole Blood", "Tissue"]
                )
            )

            st.plotly_chart(intensity_fig, use_container_width=True)

        with tab2:
            # CV plot for the filtered data
            cv_fig = px.bar(
                filtered_data_cv,
                x='Matrice',
                y='CV',
                color='Matrice',
                labels={'CV': 'Coefficient of Variation (CV)'},
                title=f"CV for {protein}"
            )
            st.plotly_chart(cv_fig, use_container_width=True)

    else:
        # Display the full data table by default or when reset is clicked
        st.subheader("Full Data Table")
        st.dataframe(data_long, use_container_width=True)  # Use full width for full data

# Run the viewer function
show_protein_total_viewer()
