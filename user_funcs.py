import numpy as np
import matplotlib.pyplot as plt
import os
import math
import Bio.PDB
import pandas as pd
import streamlit as st

def degrees(rad_angle) :
    """Converts any angle in radians to degrees.

    If the input is None, the it returns None.
    For numerical input, the output is mapped to [-180,180]
    """
    if rad_angle is None :
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180 :
        angle = angle - 360
    while angle < -180 :
        angle = angle + 360
    return angle

def ramachandran_type(residue, next_residue) :
    """Expects Bio.PDB residues, returns ramachandran 'type'

    If this is the last residue in a polypeptide, use None
    for next_residue.

    Return value is a string: "General", "Glycine", "Proline"
    or "Pre-Pro".
    """
    if residue.resname.upper()=="GLY" :
        return "Glycine"
    elif residue.resname.upper()=="PRO" :
        return "Proline"
    elif next_residue is not None \
    and next_residue.resname.upper()=="PRO" :
        #exlcudes those that are Pro or Gly
        return "Pre-Pro"
    else :
        return "General"    


def write_degrees(pdb_code):
	structure = Bio.PDB.PDBParser().get_structure(pdb_code, "%s.pdb" % pdb_code)
	output_file = open("%s_biopython.csv" % pdb_code,"w")
	for model in structure :
	    for chain in model :
	        print ("Chain %s" % str(chain.id))
	        polypeptides = Bio.PDB.CaPPBuilder().build_peptides(chain)
	        for poly_index, poly in enumerate(polypeptides) :
	            phi_psi = poly.get_phi_psi_list()
	            for res_index, residue in enumerate(poly) :
	                phi, psi = phi_psi[res_index]
	                if phi and psi :
	                    #Don't write output when missing an angle
	                    output_file.write("%s:Chain%s:%s%i,%f,%f,%s\n" \
	                        % (pdb_code, str(chain.id), residue.resname,
	                           residue.id[1], degrees(phi), degrees(psi),
	                           ramachandran_type(residue, poly[res_index+1])))
	output_file.close()


def create_figure(pdb_code):
    fig = plt.Figure()
    axis = fig.add_subplot(1,1,1)
    df = pd.read_csv(pdb_code+'_biopython.csv')
    xs = df.iloc[:,1]
    ys = df.iloc[:,2]
    axis.set_facecolor('white')
    axis.set_ylim(-180,180)
    axis.set_xlim(-180,180)
    axis.set_xticks(np.arange(-180.1,180.1,45))
    axis.set_yticks(np.arange(-180.1,180.1,45))
    axis.arrow(-180,0,360,0,color='black') 
    axis.arrow(0,-180,0,360,color='black')
    axis.scatter(xs,ys,c='r')
    st.pyplot(fig)