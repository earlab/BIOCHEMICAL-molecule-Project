/*









<PDBx:atom_siteCategory>
      <PDBx:atom_site id="1">
         <PDBx:B_iso_or_equiv>40.56</PDBx:B_iso_or_equiv> 	Isotropic atomic displacement parameter, or equivalent isotropic
               atomic displacement parameter, B~eq~, calculated from the
               anisotropic displacement parameters.
         <PDBx:Cartn_x>125.001</PDBx:Cartn_x>				The x atom-site coordinate in angstroms specified according to
               a set of orthogonal Cartesian axes related to the cell axes as
               specified by the description given in
               _atom_sites.Cartn_transform_axes.
         <PDBx:Cartn_y>11.745</PDBx:Cartn_y>				The y atom-site coordinate in angstroms specified according to
               a set of orthogonal Cartesian axes related to the cell axes as
               specified by the description given in
               _atom_sites.Cartn_transform_axes.				
         <PDBx:Cartn_z>18.819</PDBx:Cartn_z>				The z atom-site coordinate in angstroms specified according to
               a set of orthogonal Cartesian axes related to the cell axes as
               specified by the description given in
               _atom_sites.Cartn_transform_axes.		
         <PDBx:auth_asym_id>A</PDBx:auth_asym_id>			An alternative identifier for _atom_site.label_asym_id that
               may be provided by an author in order to match the identification
               used in the publication that describes the structure.					
         <PDBx:auth_atom_id>N</PDBx:auth_atom_id>			An alternative identifier for _atom_site.label_atom_id that
               may be provided by an author in order to match the identification
               used in the publication that describes the structure		
         <PDBx:auth_comp_id>ALA</PDBx:auth_comp_id>			An alternative identifier for _atom_site.label_asym_id that
               may be provided by an author in order to match the identification
               used in the publication that describes the structure.	
         <PDBx:auth_seq_id>2</PDBx:auth_seq_id>				An alternative identifier for _atom_site.label_seq_id that
               may be provided by an author in order to match the identification
               used in the publication that describes the structure.

               Note that this is not necessarily a number, that the values do
               not have to be positive, and that the value does not have to
               correspond to the value of _atom_site.label_seq_id. The value
               of _atom_site.label_seq_id is required to be a sequential list
               of positive integers.

               The author may assign values to _atom_site.auth_seq_id in any
               desired way. For instance, the values may be used to relate
               this structure to a numbering scheme in a homologous structure,
               including sequence gaps or insertion codes. Alternatively, a
               scheme may be used for a truncated polymer that maintains the
               numbering scheme of the full length polymer. In all cases, the
               scheme used here must match the scheme used in the publication
               that describes the structure.
         <PDBx:group_PDB>ATOM</PDBx:group_PDB>				The group of atoms to which the atom site belongs. This data
               item is provided for compatibility with the original Protein
               Data Bank format, and only for that purpose.
         <PDBx:label_alt_id xsi:nil="true" />				A component of the identifier for this atom site.
               For further details, see the definition of the ATOM_SITE_ALT
               category.
			
 			  Atom sites with the alternative ID set to null are not
			  modeled in alternative conformations
			  This data item is a pointer to _atom_sites_alt.id in the
               ATOM_SITES_ALT category.
         <PDBx:label_asym_id>A</PDBx:label_asym_id>			 A component of the identifier for this atom site.
               For further details, see the definition of the STRUCT_ASYM
               category.

               This data item is a pointer to _struct_asym.id in the
               STRUCT_ASYM category.
         <PDBx:label_atom_id>N</PDBx:label_atom_id>			A component of the identifier for this atom site.

               This data item is a pointer to _chem_comp_atom.atom_id in the
               CHEM_COMP_ATOM category.
         <PDBx:label_comp_id>ALA</PDBx:label_comp_id>		A component of the identifier for this atom site.

               This data item is a pointer to _chem_comp.id in the CHEM_COMP
               category.
         <PDBx:label_entity_id>1</PDBx:label_entity_id>		This data item is a pointer to _entity.id in the ENTITY category.			
         <PDBx:label_seq_id>2</PDBx:label_seq_id>			This data item is a pointer to _entity_poly_seq.num in the
               ENTITY_POLY_SEQ category.	
         <PDBx:occupancy>1.00</PDBx:occupancy>				The fraction of the atom type present at this site.
               The sum of the occupancies of all the atom types at this site
               may not significantly exceed 1.0 unless it is a dummy site.
         <PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num>	Always 1		
         <PDBx:type_symbol>N</PDBx:type_symbol>				This data item is a pointer to _atom_type.symbol in the
               ATOM_TYPE category.
      </PDBx:atom_site>
 
 		

 		Examples:

 		<PDBx:atom_siteCategory>
		      <PDBx:atom_site id="1">
		         <PDBx:B_iso_or_equiv>40.56</PDBx:B_iso_or_equiv>
		         <PDBx:Cartn_x>125.001</PDBx:Cartn_x>
		         <PDBx:Cartn_y>11.745</PDBx:Cartn_y>
		         <PDBx:Cartn_z>18.819</PDBx:Cartn_z>
		         <PDBx:auth_asym_id>A</PDBx:auth_asym_id>
		         <PDBx:auth_atom_id>N</PDBx:auth_atom_id>
		         <PDBx:auth_comp_id>ALA</PDBx:auth_comp_id>
		         <PDBx:auth_seq_id>2</PDBx:auth_seq_id>
		         <PDBx:group_PDB>ATOM</PDBx:group_PDB>
		         <PDBx:label_alt_id xsi:nil="true" />
		         <PDBx:label_asym_id>A</PDBx:label_asym_id>
		         <PDBx:label_atom_id>N</PDBx:label_atom_id>
		         <PDBx:label_comp_id>ALA</PDBx:label_comp_id>
		         <PDBx:label_entity_id>1</PDBx:label_entity_id>
		         <PDBx:label_seq_id>2</PDBx:label_seq_id>
		         <PDBx:occupancy>1.00</PDBx:occupancy>
		         <PDBx:pdbx_PDB_model_num>1</PDBx:pdbx_PDB_model_num>
		         <PDBx:type_symbol>N</PDBx:type_symbol>
		         
         
         
         
 		<PDBx:cellCategory>
		      <PDBx:cell entry_id="3SPX">
		         <PDBx:Z_PDB>4</PDBx:Z_PDB>
		         <PDBx:angle_alpha>90.000</PDBx:angle_alpha>
		         <PDBx:angle_beta>90.000</PDBx:angle_beta>
		         <PDBx:angle_gamma>90.000</PDBx:angle_gamma>
		         <PDBx:length_a>115.266</PDBx:length_a>
		         <PDBx:length_b>61.974</PDBx:length_b>
		         <PDBx:length_c>43.427</PDBx:length_c>
		      </PDBx:cell>
      
			<PDBx:chem_compCategory>
			      <PDBx:chem_comp id="ALA">
			         <PDBx:formula>C3 H7 N O2</PDBx:formula>
			         <PDBx:formula_weight>89.094</PDBx:formula_weight>
			         <PDBx:mon_nstd_flag>y</PDBx:mon_nstd_flag>
			         <PDBx:name>ALANINE</PDBx:name>
			         <PDBx:type>L-peptide linking</PDBx:type>
			      </PDBx:chem_comp>
			
<PDBx:atom_typeCategory>
      <PDBx:atom_type symbol="N"></PDBx:atom_type>
      <PDBx:atom_type symbol="C"></PDBx:atom_type>
      <PDBx:atom_type symbol="O"></PDBx:atom_type>
      <PDBx:atom_type symbol="S"></PDBx:atom_type>
      <PDBx:atom_type symbol="P"></PDBx:atom_type>
      <PDBx:atom_type symbol="CL"></PDBx:atom_type>
      <PDBx:atom_type symbol="NA"></PDBx:atom_type>			 
 
 
  </PDBx:entity_polyCategory>
   <PDBx:entity_poly_seqCategory>
      <PDBx:entity_poly_seq entity_id="1" mon_id="MET" num="1">
         <PDBx:hetero>n</PDBx:hetero>


<PDBx:pdbx_poly_seq_schemeCategory>
      <PDBx:pdbx_poly_seq_scheme asym_id="A" entity_id="1" mon_id="MET" seq_id="1">
         <PDBx:hetero>n</PDBx:hetero>
         <PDBx:ndb_seq_num>1</PDBx:ndb_seq_num>
         <PDBx:pdb_ins_code xsi:nil="true" />
         <PDBx:pdb_seq_num>1</PDBx:pdb_seq_num>
         <PDBx:pdb_strand_id>A</PDBx:pdb_strand_id>
      </PDBx:pdbx_poly_seq_scheme>


<PDBx:struct_sheet_rangeCategory>
      <PDBx:struct_sheet_range id="1" sheet_id="A">
         <PDBx:beg_auth_asym_id>A</PDBx:beg_auth_asym_id>
         <PDBx:beg_auth_comp_id>ALA</PDBx:beg_auth_comp_id>
         <PDBx:beg_auth_seq_id>24</PDBx:beg_auth_seq_id>
         <PDBx:beg_label_asym_id>A</PDBx:beg_label_asym_id>
         <PDBx:beg_label_comp_id>ALA</PDBx:beg_label_comp_id>
         <PDBx:beg_label_seq_id>24</PDBx:beg_label_seq_id>
         <PDBx:end_auth_asym_id>A</PDBx:end_auth_asym_id>
         <PDBx:end_auth_comp_id>TYR</PDBx:end_auth_comp_id>
         <PDBx:end_auth_seq_id>26</PDBx:end_auth_seq_id>
         <PDBx:end_label_asym_id>A</PDBx:end_label_asym_id>
         <PDBx:end_label_comp_id>TYR</PDBx:end_label_comp_id>
         <PDBx:end_label_seq_id>26</PDBx:end_label_seq_id>
      </PDBx:struct_sheet_range>



 <PDBx:struct_confCategory>
      <PDBx:struct_conf id="HELX_P1">
         <PDBx:beg_auth_asym_id>A</PDBx:beg_auth_asym_id>
         <PDBx:beg_auth_comp_id>SER</PDBx:beg_auth_comp_id>
         <PDBx:beg_auth_seq_id>14</PDBx:beg_auth_seq_id>
         <PDBx:beg_label_asym_id>A</PDBx:beg_label_asym_id>
         <PDBx:beg_label_comp_id>SER</PDBx:beg_label_comp_id>
         <PDBx:beg_label_seq_id>14</PDBx:beg_label_seq_id>
         <PDBx:conf_type_id>HELX_P</PDBx:conf_type_id>
         <PDBx:end_auth_asym_id>A</PDBx:end_auth_asym_id>
         <PDBx:end_auth_comp_id>ILE</PDBx:end_auth_comp_id>
         <PDBx:end_auth_seq_id>19</PDBx:end_auth_seq_id>
         <PDBx:end_label_asym_id>A</PDBx:end_label_asym_id>
         <PDBx:end_label_comp_id>ILE</PDBx:end_label_comp_id>
         <PDBx:end_label_seq_id>19</PDBx:end_label_seq_id>
         <PDBx:pdbx_PDB_helix_class>5</PDBx:pdbx_PDB_helix_class>
         <PDBx:pdbx_PDB_helix_id>1</PDBx:pdbx_PDB_helix_id>
         <PDBx:pdbx_PDB_helix_length>6</PDBx:pdbx_PDB_helix_length>
      </PDBx:struct_conf>


*/
/*
Alanine
Ala
A
nonpolar
neutral
1.8


Arginine
Arg
R
polar
positive
−4.5


Asparagine
Asn
N
polar
neutral
−3.5


Aspartic acid
Asp
D
polar
negative
−3.5


Cysteine
Cys
C
polar
neutral
2.5
250
0.3
Glutamic acid
Glu
E
polar
negative
−3.5


Glutamine
Gln
Q
polar
neutral
−3.5


Glycine
Gly
G
nonpolar
neutral
−0.4


Histidine
His
H
polar
positive(10%)
neutral(90%)
−3.2
211
5.9
Isoleucine
Ile
I
nonpolar
neutral
4.5


Leucine
Leu
L
nonpolar
neutral
3.8


Lysine
Lys
K
polar
positive
−3.9


Methionine
Met
M
nonpolar
neutral
1.9


Phenylalanine
Phe
F
nonpolar
neutral
2.8
257, 206, 188
0.2, 9.3, 60.0
Proline
Pro
P
nonpolar
neutral
−1.6


Serine
Ser
S
polar
neutral
−0.8


Threonine
Thr
T
polar
neutral
−0.7


Tryptophan
Trp
W
nonpolar
neutral
−0.9
280, 219
5.6, 47.0
Tyrosine
Tyr
Y
polar
neutral
−1.3
274, 222, 193
1.4, 8.0, 48.0
Valine
Val
V
nonpolar
neutral
4.2





*/

