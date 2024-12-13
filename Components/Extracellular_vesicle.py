import streamlit as st
import pandas as pd
import os
import plotly.express as px

def show_protein_extracellular_vesicle():
    # --- Data Loading and Processing ---
    
    # Load the new data file
    data_path = os.path.join('data', 'EV_protein_list.csv')
    data = pd.read_csv(data_path)

    # Truncate column names by removing ".d.PG.Quantity"
    data.columns = data.columns.str.replace(r'\\.d\\.PG\\.Quantity$', '', regex=True)

    # --- Streamlit App Interface ---
    st.title("Extracellular Vesicles")

    # Sidebar inputs for protein search and selection
    st.sidebar.header("Search Options")
    protein_options = data['PG.ProteinAccessions'].unique()
    protein_search = st.sidebar.text_input("Enter PG.ProteinAccessions", key="filtered_protein_search")
    protein_select = st.sidebar.multiselect("Select PG.ProteinAccessions", protein_options, key="filtered_protein_select")
    search = st.sidebar.button("Search", key="filtered_search_button")
    reset = st.sidebar.button("Reset to Full Table", key="filtered_reset_button")

    # Display either full or filtered data based on search or reset
    if search:
        # Combine entered proteins and selected proteins
        entered_proteins = [p.strip() for p in protein_search.split(',') if p.strip()]
        selected_proteins = protein_select
        proteins_to_filter = list(set(entered_proteins + selected_proteins))

        filtered_data = data[data['PG.ProteinAccessions'].isin(proteins_to_filter)]
        st.subheader(f"Data Table for {', '.join(proteins_to_filter)}")
        st.dataframe(filtered_data, use_container_width=True)  # Use full width for filtered data

        if len(proteins_to_filter) > 5:
            st.warning("The bar chart currently supports displaying data for up to 5 proteins. Please reduce the number of selected proteins.")
        else:
            # Prepare data for plotting
            melted_data = filtered_data.melt(
                id_vars=['PG.ProteinAccessions', 'PG.Genes', 'PG.ProteinDescriptions', 'PG.ProteinNames', 'PG.OrganismId'],
                var_name='Sample',
                value_name='Intensity'
            )

            # Display the MS Quantity Plot
            st.subheader("MS Intensity Plot")
            quantity_fig = px.bar(
                melted_data,
                x='Sample',
                y='Intensity',
                color='PG.ProteinAccessions',
                labels={'Intensity': 'MS Intensity', 'Sample': 'Sample'},
                title="MS Intensity for Selected Proteins",
                barmode='group'
            )

            # Update layout to align categories
            quantity_fig.update_layout(
                xaxis=dict(
                    title='Sample',
                    categoryorder="array",
                    categoryarray=[
                        '[1] EV_071A_R1_S2-A1_1_9412',
                        '[2] EV_071A_R2_S2-A2_1_9413',
                        '[3] EV_071A_R3_S2-A3_1_9414',
                        '[4] EV_071A_R4_S2-A4_1_9415'
                    ]
                )
            )

            st.plotly_chart(quantity_fig, use_container_width=True)

    else:
        # Display the full data table by default or when reset is clicked
        st.subheader("Full Data Table")
        st.dataframe(data, use_container_width=True)  # Use full width for full data

# Run the viewer function
show_protein_extracellular_vesicle()
