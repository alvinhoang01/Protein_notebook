import streamlit as st
import pandas as pd
import numpy as np
import os
import plotly.express as px

def show_protein_viewer():
    # --- Data Loading and Processing ---

    # Read the data from the 'data' folder (adjust the path accordingly)
    data_path = os.path.join('data', 'method_total.csv')
    data = pd.read_csv(data_path)

    # Pivot the data into long format
    data_long = pd.melt(data, id_vars=['Protein'], var_name='Sample', value_name='Intensity')

    # Create 'Method' and 'Replicate' columns based on patterns in 'key'
    data_long['Method'] = np.select(
        [
            data_long['Sample'].str.contains('Preomics'),
            data_long['Sample'].str.contains('Thermokit'),
            data_long['Sample'].str.contains('Direct'),
            data_long['Sample'].str.contains('Seer'),
            data_long['Sample'].str.contains('CO'),
            data_long['Sample'].str.contains('EV'),
            data_long['Sample'].str.contains('SPEG')
        ],
        ['Preomics', 'Thermo', 'Direct', 'Seer', 'CO', 'EV', 'SPEG'],
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

    # Map Replicate to specific colors for the legend
    replicate_color_map = {'R1': 'red', 'R2': 'blue', 'R3': 'green', 'R4': 'purple'}
    data_long['ReplicateColor'] = data_long['Replicate'].map(replicate_color_map)

    # Calculate mean, standard deviation, and CV for each Protein-Method group
    data_cv = data_long.groupby(['Protein', 'Method']).agg(
        mean_Intensity=('Intensity', 'mean'),
        sd_Intensity=('Intensity', 'std')
    ).reset_index()
    data_cv['CV'] = (data_cv['sd_Intensity'] / data_cv['mean_Intensity']) * 100

    # Merge the CV back into the long data frame
    data_long = pd.merge(data_long, data_cv[['Protein', 'Method', 'mean_Intensity', 'sd_Intensity', 'CV']], on=['Protein', 'Method'], how='left')

    # --- Streamlit App Interface ---
    st.title("Methods")

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
        st.dataframe(filtered_data, use_container_width=True)  # Always display the data table

        if len(proteins_to_filter) > 4:
            st.warning("The current data visualization supports only up to 4 proteins. Please reduce your selection for visualization.")
        else:
            filtered_data_cv = data_cv[data_cv['Protein'].isin(proteins_to_filter)]

            # Create tabs for the two plots
            tab1, tab2 = st.tabs(["MS Intensity Plot", "CV Plot"])

            with tab1:
                # Create side-by-side bar charts for each selected protein
                col1, col2 = st.columns(2)
                for i, protein in enumerate(proteins_to_filter):
                    protein_data = filtered_data[filtered_data['Protein'] == protein]
                    mean_data = protein_data.drop_duplicates(subset=['Method', 'mean_Intensity'])

                    with (col1 if i % 2 == 0 else col2):
                        # st.markdown(f"### {protein}")
                        intensity_fig = px.bar(
                            mean_data,
                            x='Method',
                            y='mean_Intensity',
                            labels={'mean_Intensity': 'Average MS Intensity'},
                            title=f"MS Intensity for {protein}"
                        )
                        for replicate, color in replicate_color_map.items():
                            replicate_data = protein_data[protein_data['Replicate'] == replicate]
                            intensity_fig.add_scatter(
                                x=replicate_data['Method'],
                                y=replicate_data['Intensity'],
                                mode='markers',
                                marker=dict(size=8, color=color),
                                name=replicate
                            )
                        intensity_fig.update_layout(legend_title_text="Replicates")
                        st.plotly_chart(intensity_fig, use_container_width=True)

            with tab2:
                # Ensure methods in the filtered_data_cv are ordered correctly
                method_order = ['Preomics', 'Thermo', 'Direct', 'Seer', 'CO', 'EV', 'SPEG']
                filtered_data_cv['Method'] = pd.Categorical(filtered_data_cv['Method'], categories=method_order, ordered=True)

                col1, col2 = st.columns(2)
                for i, protein in enumerate(proteins_to_filter):
                    protein_data_cv = filtered_data_cv[filtered_data_cv['Protein'] == protein]

                    # Sort by Method to ensure proper ordering
                    protein_data_cv = protein_data_cv.sort_values(by='Method')

                    with (col1 if i % 2 == 0 else col2):
                        # Create the CV bar plot
                        cv_fig = px.bar(
                            protein_data_cv,
                            x='Method',
                            y='CV',
                            labels={'CV': 'Coefficient of Variation (CV)'},
                            title=f"CV for {protein}"
                        )
                        cv_fig.update_layout(legend_title_text="Replicates")
                        st.plotly_chart(cv_fig, use_container_width=True)

    else:
        # Display the full data table by default or when reset is clicked
        st.subheader("Full Data Table")
        st.dataframe(data_long, use_container_width=True)  # Use full width for full data

# Run the viewer function
show_protein_viewer()
