q.put(  // Slower Method x10 Slower than torsionAngles
~torAngles,{|env, seqId|
	var carCAp, carCp, carN, carCA, carC, carNa;				//Cartesian
	var phi, psi, omega;
	var nodes;
	nodes = Array.newClear(6);
	nodes.put(0,
		env.pdbml.getDocumentElement.getElement("PDBx:atom_siteCategory").getChildNodes
		.detect {|item|
		item.getElement("PDBx:label_seq_id").getText.asInteger == (seqId-1) &&
		item.getElement("PDBx:label_atom_id").getText == "CA"
		};
	);
	nodes.put(1, nodes.at(0).getNextSibling);
	nodes.put(2,
		env.pdbml.getDocumentElement.getElement("PDBx:atom_siteCategory").getChildNodes
		.detect {|item|
		item.getElement("PDBx:label_seq_id").getText.asInteger == (seqId) &&
		item.getElement("PDBx:label_atom_id").getText == "N"
		};
	);
	nodes.put(3, nodes.at(2).getNextSibling);
	nodes.put(4, nodes.at(3).getNextSibling);
	nodes.put(5,
		env.pdbml.getDocumentElement.getElement("PDBx:atom_siteCategory").getChildNodes
		.detect {|item|
		item.getElement("PDBx:label_seq_id").getText.asInteger == (seqId+1) &&
		item.getElement("PDBx:label_atom_id").getText == "N"
		};
	);
carCAp = Cartesian(
		x: nodes[0].getElement("PDBx:Cartn_x").getText.asFloat,
		y: nodes[0].getElement("PDBx:Cartn_y").getText.asFloat,
		z: nodes[0].getElement("PDBx:Cartn_z").getText.asFloat);

carCp = Cartesian(
		x: nodes[1].getElement("PDBx:Cartn_x").getText.asFloat,
		y: nodes[1].getElement("PDBx:Cartn_y").getText.asFloat,
		z: nodes[1].getElement("PDBx:Cartn_z").getText.asFloat);

carN = Cartesian(
		x: nodes[2].getElement("PDBx:Cartn_x").getText.asFloat,
		y: nodes[2].getElement("PDBx:Cartn_y").getText.asFloat,
		z: nodes[2].getElement("PDBx:Cartn_z").getText.asFloat);

carCA = Cartesian(
		x: nodes[3].getElement("PDBx:Cartn_x").getText.asFloat,
		y: nodes[3].getElement("PDBx:Cartn_y").getText.asFloat,
		z: nodes[3].getElement("PDBx:Cartn_z").getText.asFloat);

carC = Cartesian(
		x: nodes[4].getElement("PDBx:Cartn_x").getText.asFloat,
		y: nodes[4].getElement("PDBx:Cartn_y").getText.asFloat,
		z: nodes[4].getElement("PDBx:Cartn_z").getText.asFloat);

carNa = Cartesian(
		x: nodes[5].getElement("PDBx:Cartn_x").getText.asFloat,
		y: nodes[5].getElement("PDBx:Cartn_y").getText.asFloat,
		z: nodes[5].getElement("PDBx:Cartn_z").getText.asFloat);

phi = env.torFromCart(carCp, carN, carCA, carC);
psi = env.torFromCart(carN, carCA, carC, carNa);
omega = env.torFromCart(carCAp, carCp, carN, carCA);
[phi, psi, omega];
};
)