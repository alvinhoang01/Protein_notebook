import streamlit as st
import os

def main():
    st.markdown(
        """
        <style>
            .title {
                font-family: "Arial", sans-serif;
                color: #008080;
                font-size: 42px;
                font-weight: bold;
            }
            .subtitle {
                font-family: "Arial", sans-serif;
                color: #800000;
                font-size: 22px;
            }
            .main-header {
                font-family: "Arial", sans-serif; 
                font-size: 28px;
                font-weight: bold;
                text-decoration: underline;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("<div class='title'>Human Blood Proteomics Atlas</div>", unsafe_allow_html=True)
    st.write("\n")

    # Display the static image (logo or decorative image)
    data_path = os.path.join('data', "Protein_logo.png")
    st.image(data_path, use_column_width=True)

    st.write("""
    Our web platform enables researchers to access mass spectrometry-based proteomics data on protein detectability, signal intensity, and quantification reproducibility across various sample preparation methods. Explore features including:
    - Comparative performance analysis of blood proteins across different sample preparation techniques
    - Comprehensive performance analysis of proteins in different clinical samples
    """)

    # Example of an exit button
    left_column2, right_column2 = st.columns([18, 2])
    with right_column2:
        st.write("\n")
        st.write("\n")
        if st.button('Exit'):
            st.error("The app has been stopped. Please close the browser window manually.")
            os._exit(0)

# This file will be imported by `app.py` to render the home page
