import streamlit as st
import os

def main():
    st.markdown(
        """
        <style>
            .title {
                font-family: "Arial", sans-serif;
                color: #008080;
                font-size: 36px;
                font-weight: bold;
            }
            .subtitle {
                font-family: "Arial", sans-serif;
                color: #800000;
                font-size: 35px;
            }
            .main-header {
                font-family: "Arial", sans-serif; 
                font-size: 35px;
                font-weight: bold;
                text-decoration: underline;
            }
            .description {
                font-family: "Arial", sans-serif;
                font-size: 22px;  
                line-height: 1.6;
            }
            .description li {
                font-size: 19px;
                margin-bottom: 8px;
            }
        </style>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("<div class='title'>Human Blood Proteomics Atlas</div>", unsafe_allow_html=True)
    st.write("\n")

    # Display the static image (logo or decorative image)
    from PIL import Image

    data_path = os.path.join('data', "Workflow_V4.tif")
    image = Image.open(data_path)

    # Create three columns and put the image in the center one
    col1, col2, col3 = st.columns([1, 2, 2])
    with col2:
        st.image(image, width=800)


    st.markdown(
    """
    <div class='description'>
        Our web platform enables researchers to access mass spectrometry-based proteomics data on protein detectability, signal intensity, and quantification reproducibility across various sample preparation methods. Explore features including:
        <ul>
            <li>Comparative performance analysis of blood proteins across different sample preparation techniques</li>
            <li>Comprehensive performance analysis of proteins in different clinical samples</li>
        </ul>
    </div>
    """,
    unsafe_allow_html=True
)

    # Example of an exit button
    
    left_column2, right_column2 = st.columns([18, 2])
    with right_column2:
        st.write("\n")
        st.write("\n")
        if st.button('Exit'):
            st.error("The app has been stopped. Please close the browser window manually.")
            os._exit(0)
    

# This file will be imported by `app.py` to render the home page
