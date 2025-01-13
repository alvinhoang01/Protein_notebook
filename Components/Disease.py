import streamlit as st
import pandas as pd
import os
import plotly.express as px

def show_protein_disease():
    # --- Data Loading and Processing ---
    
    # Load the new data file
    data_path = os.path.join('data', 'Normal_disease_protein_total.csv')
    data = pd.read_csv(data_path)

    # Replace values for better readability
    data = data.replace(
        {1: 'Present', 2: 'Present', 3: 'Present', 'N/A': 'Absent'}
    )

    # --- Streamlit App Interface ---
    st.title("Disease")

    # Sidebar inputs for protein search and selection
    st.sidebar.header("Search Options")
    protein_options = data['PG.ProteinAccessions'].unique()
    protein_search = st.sidebar.text_input("Enter PG.ProteinAccessions", key="protein_search")
    protein_select = st.sidebar.multiselect("Select PG.ProteinAccessions", protein_options, key="protein_select")
    search = st.sidebar.button("Search", key="search_button")
    reset = st.sidebar.button("Reset", key="reset_button")

    # Filter data based on search or selection
    if search:
        entered_proteins = [p.strip() for p in protein_search.split(',') if p.strip()]
        selected_proteins = protein_select
        proteins_to_filter = list(set(entered_proteins + selected_proteins))

        filtered_data = data[data['PG.ProteinAccessions'].isin(proteins_to_filter)]
        st.subheader(f"Filtered Data for {', '.join(proteins_to_filter)}")
        st.dataframe(filtered_data, use_container_width=True)

        # Bar chart summarizing presence across sample types
        presence_summary = filtered_data.melt(id_vars=['PG.ProteinAccessions'], var_name='Sample_Type', value_name='Presence')
        presence_summary = presence_summary.groupby(['Sample_Type', 'Presence']).size().reset_index(name='Count')
        bar_chart = px.bar(
            presence_summary, 
            x='Sample_Type', 
            y='Count', 
            color='Presence', 
            barmode='group',
            title="Presence Summary Across Sample Types"
        )
        bar_chart.update_traces(width=0.4)  # Make bars thinner
        bar_chart.update_layout(showlegend=False)  # Remove legend
        st.plotly_chart(bar_chart, use_container_width=True)
    else:
        # Display full data by default or on reset
        st.subheader("Full Data Table (Presence/Absence)")
        st.dataframe(data, use_container_width=True)

# Run the viewer function
show_protein_disease()
