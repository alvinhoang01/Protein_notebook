import streamlit as st
from datetime import datetime
from streamlit_option_menu import option_menu

# --- Display Sidebar and Pages ---
def display_sidebar_and_pages():
    st.sidebar.title("Welcome!")
    st.sidebar.markdown("#### Select a page to navigate:")

    with st.sidebar:
        page = option_menu(
            menu_title=None,
            options=["Home Page", "Protein Notebook", "Matrix Analysis"],
            icons=["house", "database", "bar-chart-line"],
            menu_icon="cast",
            default_index=0,
            orientation="vertical",
            styles={
                "nav-link": {
                    "font-size": "16px",
                    "text-align": "left",
                    "margin": "5px",
                    "--hover-color": "#f0f0f0",
                },
                "icon": {
                    "font-size": "18px",
                    "color": "#008cba",
                },
                "container": {
                    "padding": "5px",
                    "background-color": "transparent"
                },
                "nav-link-selected": {
                    "background-color": "rgba(0, 123, 255, 0.15)",
                    "font-weight": "bold",
                    "color": "#008cba",
                },
            },
        )

    # --- Load Pages Based on Sidebar Selection ---
    if page == "Home Page":
        from Home_page import main as show_home_page
        show_home_page()

    elif page == "Protein Notebook":
        from ptmdatabase.Components.protein_notebook import show_protein_viewer as show_database_page
        show_database_page()

    elif page == "Matrix Analysis":
        from Components.Matrix_analysis import main as show_matrix_analysis_page
        show_matrix_analysis_page()

    # No logout button needed since authentication is removed
    # if st.sidebar.button('🔒 Logout'):
    #     st.session_state['authenticated'] = False
    #     st.session_state['username'] = ''
    #     st.session_state['signup'] = False

# --- Main Function ---
def main():
    st.set_page_config(
        page_title="Proteoform Database Generation App",
        page_icon="🧬",
        layout="wide",
        initial_sidebar_state="collapsed",
    )

    # Directly display the main content without authentication
    display_sidebar_and_pages()

if __name__ == '__main__':
    main()
