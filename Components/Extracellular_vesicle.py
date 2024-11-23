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
    data.columns = data.columns.str.replace(r'\.d\.PG\.Quantity$', '', regex=True)

    # --- Streamlit App Interface ---
    st.title("Extracellular vesicles")

    # Sidebar inputs for protein search and selection
    st.sidebar.header("Search Options")
    protein_options = data['PG.ProteinAccessions'].unique()
    protein_search = st.sidebar.text_input("Enter PG.ProteinAccessions", key="filtered_protein_search")
    protein_select = st.sidebar.selectbox("Select PG.ProteinAccessions", protein_options, key="filtered_protein_select")
    search = st.sidebar.button("Search", key="filtered_search_button")
    reset = st.sidebar.button("Reset to Full Table", key="filtered_reset_button")

    # Display either full or filtered data based on search or reset
    if search:
        # Determine which protein to filter by
        protein = protein_search if protein_search != "" else protein_select
        filtered_data = data[data['PG.ProteinAccessions'] == protein]

        # Display filtered data and plots
        st.subheader(f"Data Table for {protein}")
        st.dataframe(filtered_data, use_container_width=True)  # Use full width for filtered data

        # Display the MS Quantity Plot
        st.subheader("MS Intensity Plot")
        quantity_fig = px.bar(
            filtered_data.melt(id_vars=['PG.ProteinAccessions', 'PG.Genes', 'PG.ProteinDescriptions', 'PG.ProteinNames', 'PG.OrganismId']),
            x='variable',  # Melted columns as x-axis
            y='value',     # Quantity values as y-axis
            color='variable',
            labels={'value': 'Intensity', 'variable': 'Sample'},
            title=f"MS Intensity for {protein}"
        )
        
        # Update layout to remove legend and align categories
        quantity_fig.update_layout(
            showlegend=True,  # Hide legend
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
