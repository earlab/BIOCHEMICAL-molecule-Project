d = DOMDocument("/Users/stratosbichakis/Downloads/components-pub-xml/ARG.xml");
//categories listing
d.getDocumentElement.getChildNodes.do({|item| item.getNodeName.postln;})

~chem_comp_atomCat = d.getDocumentElement.getElement("PDBx:chem_comp_atomCategory");
~chem_comp_bondCat = d.getDocumentElement.getElement("PDBx:chem_comp_bondCategory");

~chem_comp_atomCat.getFirstChild.getAttributes;
~chem_comp_bondCat.getFirstChild.getAttributes;

//Build chem_comp_atom Dictionary
~chem_comp_atomDic = ~chem_comp_atomCat.getChildNodes.collectAs({|item| item.getAttribute("atom_id").asSymbol -> item.getChildNodes.collectAs({|item| item.getNodeName.asSymbol-> item.getText}, Dictionary)}, Dictionary);

//Build chem_comp_bond List of connected pairs
~chem_comp_bondList = ~chem_comp_bondCat.getChildNodes.collectAs({|item| [item.getAttribute("atom_id_1").asSymbol,item.getAttribute("atom_id_2").asSymbol]}, List)

//Build chem_comp_bond Dictionary
~chem_comp_bondDict = ~chem_comp_bondCat.getChildNodes.collectAs({|item| item.getChildNodes.collectAs({|item| item.getNodeName.asSymbol->item.getText}, Dictionary)}, List)

//Calculation Functions
(
~angleFromCart = {|a, b, c|
	acos(((a-b)*(c-b)).asArray.sum/(dist(a,b)*dist(b,c)))*180/pi
};

~crossFromCart = {|a, b, c|
	var e, f;
		e = (b-a).asArray *.t (c-b).asArray;
		f = [e[1][2]-e[2][1], e[2][0]-e[0][2], e[0][1]-e[1][0]];
};

~torFromCart = {|a,b,c,d|
	var e,f,g,t;
	e = ~crossFromCart.(a, b, c);
	f = ~crossFromCart.(d, c, b);
	g = (e*f).sum/(sqrt(e.collect(_.squared).sum)*sqrt(f.collect(_.squared).sum));
	t = acos(g)*180/pi;
};
)


(
a = Cartesian(
	x: ~chem_comp_atomDic['N']['PDBx:model_Cartn_x'].asFloat,
	y: ~chem_comp_atomDic['N']['PDBx:model_Cartn_y'].asFloat,
	z: ~chem_comp_atomDic['N']['PDBx:model_Cartn_z'].asFloat);

b = Cartesian(
	x: ~chem_comp_atomDic['CA']['PDBx:model_Cartn_x'].asFloat,
	y: ~chem_comp_atomDic['CA']['PDBx:model_Cartn_y'].asFloat,
	z: ~chem_comp_atomDic['CA']['PDBx:model_Cartn_z'].asFloat);

c = Cartesian(
	x: ~chem_comp_atomDic['CB']['PDBx:model_Cartn_x'].asFloat,
	y: ~chem_comp_atomDic['CB']['PDBx:model_Cartn_y'].asFloat,
	z: ~chem_comp_atomDic['CB']['PDBx:model_Cartn_z'].asFloat);

d = Cartesian(
	x: ~chem_comp_atomDic['CG']['PDBx:model_Cartn_x'].asFloat,
	y: ~chem_comp_atomDic['CG']['PDBx:model_Cartn_y'].asFloat,
	z: ~chem_comp_atomDic['CG']['PDBx:model_Cartn_z'].asFloat);
)

//Calculate Distance
dist(a, b).round(0.001);

//Calculate the bond angle
~angleFromCart.(a, b, c).round(0.1);

//Calculate torsion angle
~torFromCart.(a, b, c, d).round(0.1)+180;
)

//Number of atoms
~chem_comp_atomDic.size

//Chiral Atoms
~chem_comp_atomDic.keys.select({|item| ~chem_comp_atomDic[item]['PDBx:pdbx_stereo_config']=="S"})

//Leaving Atoms
~chem_comp_atomDic.keys.select({|item| ~chem_comp_atomDic[item]['PDBx:pdbx_leaving_atom_flag']=="Y"})

//Number of bonds
~chem_comp_bondList.size

//Test if bond is valid
~chem_comp_bondList.any({|item| item.includesAll(['CA','CB'])})
~chem_comp_bondList.detectIndex({|item| item.includesAll(['CA','CB'])})

//Test if there is any aromatic bont
~chem_comp_bondDict.any({|item| item['PDBx:pdbx_aromatic_flag']=="Y"})