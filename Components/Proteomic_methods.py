import streamlit as st
import pandas as pd
import numpy as np
import os
import plotly.express as px

def show_protein_viewer():
    # --- Data Loading and Processing ---
    
    # Read the data from the 'data' folder (adjust the path accordingly)
    data_path = os.path.join('data', 'Protein_filter.csv')
    data = pd.read_csv(data_path)

    # Pivot the data into long format
    data_long = pd.melt(data, id_vars=['Protein'], var_name='key', value_name='value')

    # Create 'Method' and 'Replicate' columns based on patterns in 'key'
    data_long['Method'] = np.select(
        [
            data_long['key'].str.contains('Preomics'),
            data_long['key'].str.contains('Thermokit'),
            data_long['key'].str.contains('Direct'),
            data_long['key'].str.contains('Seer_NPA'),
            data_long['key'].str.contains('Seer_NPB')
        ],
        ['Preomics', 'Thermo', 'Direct', 'Seer_NPA', 'Seer_NPB'],
        default=np.nan
    )

    data_long['Replicate'] = np.select(
        [
            data_long['key'].str.contains('R1'),
            data_long['key'].str.contains('R2'),
            data_long['key'].str.contains('R3'),
            data_long['key'].str.contains('R4')
        ],
        ['R1', 'R2', 'R3', 'R4'],
        default=np.nan
    )

    # Calculate mean, standard deviation, and CV for each Protein-Method group
    data_cv = data_long.groupby(['Protein', 'Method']).agg(
        mean_value=('value', 'mean'),
        sd_value=('value', 'std')
    ).reset_index()
    data_cv['CV'] = (data_cv['sd_value'] / data_cv['mean_value']) * 100

    # Merge the CV back into the long data frame
    data_long = pd.merge(data_long, data_cv[['Protein', 'Method', 'CV']], on=['Protein', 'Method'], how='left')

    # --- Streamlit App Interface ---
    st.title("Proteomic Methods")

    # Sidebar inputs for protein search and selection
    st.sidebar.header("Search Options")
    protein_search = st.sidebar.text_input("Enter PG.ProteinAccessions", key="protein_search_method")
    protein_options = data_long['Protein'].unique()
    protein_select = st.sidebar.selectbox("Select PG.ProteinAccessions", protein_options)
    search = st.sidebar.button("Search")
    reset = st.sidebar.button("Reset to Full Table")

    # Display either full or filtered data based on search or reset
    if search:
        # Determine which protein to filter by
        protein = protein_search if protein_search != "" else protein_select
        filtered_data = data_long[data_long['Protein'] == protein]
        filtered_data_cv = data_cv[data_cv['Protein'] == protein]

        # If there are multiple values per method, consider aggregating before plotting
        aggregated_data = filtered_data.groupby(['Method']).agg(
            mean_intensity=('value', 'mean')
        ).reset_index()

        # Display filtered data and plots
        st.subheader(f"Data Table for {protein}")
        st.dataframe(filtered_data, use_container_width=True)  # Use full width for filtered data

        # Create tabs for the two plots
        tab1, tab2 = st.tabs(["MS Intensity Plot", "CV Plot"])

        with tab1:
            # Use aggregated data if needed
            intensity_fig = px.bar(
                aggregated_data,
                x='Method',
                y='mean_intensity',
                color='Method',
                labels={'mean_intensity': 'Mean MS Intensity'},
                title=f"MS Intensity for {protein}"
            )
            st.plotly_chart(intensity_fig, use_container_width=True)

        with tab2:
            cv_fig = px.bar(
                filtered_data_cv,
                x='Method',
                y='CV',
                color='Method',
                labels={'CV': 'Coefficient of Variation (CV)'},
                title=f"CV for {protein}"
            )
            st.plotly_chart(cv_fig, use_container_width=True)

    else:
        # Display the full data table by default or when reset is clicked
        st.subheader("Full Data Table")
        st.dataframe(data_long, use_container_width=True)  # Use full width for full data

# Run the viewer function
show_protein_viewer()
