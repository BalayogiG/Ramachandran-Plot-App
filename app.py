
import io
import os
import math
import Bio.PDB
import requests
import matplotlib.pyplot as plt
import matplotlib.patches as pat
import pandas as pd
import numpy as np
import sys
import warnings
import streamlit as st 
import user_funcs as uf
warnings.filterwarnings('ignore')

MASTER_URL = "https://files.rcsb.org/download/"
FORMAT = ".pdb"

st.title('Ramachandram Plot Application')
st.subheader('Build with Streamlit and BioPython')

pdb_code = st.sidebar.text_input('Enter PDB code:')

if st.button('Show Plot'):
	pdb_code = pdb_code.upper()
	# Getting data from URL
	url = MASTER_URL + pdb_code + FORMAT
	r = requests.get(url)
	with open("file.pdb", 'wb') as f:
		f.write(r.content)
	if pdb_code == "":
		st.write("Enter PDB code and press Enter....")
	else:
		uf.write_degrees("file")
		uf.create_figure('file')

