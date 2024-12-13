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
    data_long = pd.melt(data, id_vars=['Protein'], var_name='Sample', value_name='Intensity')

    # Create 'Method' and 'Replicate' columns based on patterns in 'key'
    data_long['Method'] = np.select(
        [
            data_long['Sample'].str.contains('Preomics'),
            data_long['Sample'].str.contains('Thermokit'),
            data_long['Sample'].str.contains('Direct'),
            data_long['Sample'].str.contains('Seer_NPA'),
            data_long['Sample'].str.contains('Seer_NPB')
        ],
        ['Preomics', 'Thermo', 'Direct', 'Seer_NPA', 'Seer_NPB'],
        default=np.nan
    )

    data_long['Replicate'] = np.select(
        [
            data_long['Sample'].str.contains('R1'),
            data_long['Sample'].str.contains('R2'),
            data_long['Sample'].str.contains('R3'),
            data_long['Sample'].str.contains('R4')
        ],
        ['R1', 'R2', 'R3', 'R4'],
        default=np.nan
    )

    # Calculate mean, standard deviation, and CV for each Protein-Method group
    data_cv = data_long.groupby(['Protein', 'Method']).agg(
        mean_Intensity=('Intensity', 'mean'),
        sd_Intensity=('Intensity', 'std')
    ).reset_index()
    data_cv['CV'] = (data_cv['sd_Intensity'] / data_cv['mean_Intensity']) * 100

    # Merge the CV back into the long data frame
    data_long = pd.merge(data_long, data_cv[['Protein', 'Method', 'CV']], on=['Protein', 'Method'], how='left')

    # --- Streamlit App Interface ---
    st.title("Proteomic Methods")

    # Sidebar inputs for protein search and selection
    st.sidebar.header("Search Options")
    protein_options = data_long['Protein'].unique()
    protein_search = st.sidebar.text_input("Enter PG.ProteinAccessions", key="total_protein_search")
    protein_select = st.sidebar.multiselect("Select PG.ProteinAccessions", protein_options, key="total_protein_select")
    search = st.sidebar.button("Search", key="total_search_button")
    reset = st.sidebar.button("Reset to Full Table", key="total_reset_button")

    # Display either full or filtered data based on search or reset
    if search:
        # Combine entered proteins and selected proteins
        entered_proteins = [p.strip() for p in protein_search.split(',') if p.strip()]
        selected_proteins = protein_select
        proteins_to_filter = list(set(entered_proteins + selected_proteins))

        filtered_data = data_long[data_long['Protein'].isin(proteins_to_filter)]
        st.subheader(f"Data Table for {', '.join(proteins_to_filter)}")
        st.dataframe(filtered_data, use_container_width=True)  # Use full width for filtered data

        if len(proteins_to_filter) > 5:
            st.warning("The bar chart currently supports displaying data for up to 5 proteins. Please reduce the number of selected proteins.")
        else:
            filtered_data_cv = data_cv[data_cv['Protein'].isin(proteins_to_filter)]

            # If there are multiple values per method, consider aggregating before plotting
            aggregated_data = filtered_data.groupby(['Method', 'Protein']).agg(
                mean_intensity=('Intensity', 'mean')
            ).reset_index()

            # Create tabs for the two plots
            tab1, tab2 = st.tabs(["MS Intensity Plot", "CV Plot"])

            with tab1:
                intensity_fig = px.bar(
                    aggregated_data,
                    x='Method',
                    y='mean_intensity',
                    color='Protein',
                    barmode='group',
                    labels={'mean_intensity': 'Mean MS Intensity'},
                    title=f"MS Intensity for Selected Proteins"
                )
                st.plotly_chart(intensity_fig, use_container_width=True)

            with tab2:
                cv_fig = px.bar(
                    filtered_data_cv,
                    x='Method',
                    y='CV',
                    color='Protein',
                    barmode='group',
                    labels={'CV': 'Coefficient of Variation (CV)'},
                    title=f"CV for Selected Proteins"
                )
                st.plotly_chart(cv_fig, use_container_width=True)

    else:
        # Display the full data table by default or when reset is clicked
        st.subheader("Full Data Table")
        st.dataframe(data_long, use_container_width=True)  # Use full width for full data

# Run the viewer function
show_protein_viewer()