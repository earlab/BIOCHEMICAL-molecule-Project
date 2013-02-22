// Getting documentation from PDBML Schema

h = DOMDocument("/Users/stratosbichakis/Downloads/pdbx-v40-v4.025.xsd");

//Get All the Category names
h.getDocumentElement.getChildNodes.do({|item| item.getAttribute("name").postln;})

//Get Documentation on the Category
h.getDocumentElement.getChildNodes.detect({|item| item.getAttribute("name")== ("chem_comp_atom" ++ "Type")}).getElement("xsd:annotation").getLastChild.getText

//or
h.getDocumentElement.getChildNodes.detect({|item| item.getAttribute("name")== "database_PDB_matrixType"}).getElement("xsd:annotation").getElement("xsd:documentation").getText

//Get more specific documentation
h.getDocumentElement.getChildNodes.detect({|item| item.getAttribute("name")== "atom_siteType"}).getElement("xsd:sequence").getElement("xsd:element").getElement("xsd:complexType").getElement("xsd:all").getChildNodes.detect({|item| item.getAttribute("name")=="aniso_U12_esd"}).getElement("xsd:annotation").getElement("xsd:documentation").getText

// or
h.getElementsByTagName("xsd:element").detect({|item| item.getAttribute("name")=="entity_poly_seq"}).getElement("xsd:annotation").getElement("xsd:documentation").getText


(
h = DOMDocument("/Users/stratosbichakis/Downloads/pdbx-v40-v4.025.xsd");
~help = {arg what;
h.getElementsByTagName("xsd:element").detect({|item| item.getAttribute("name")==what}).getElement("xsd:annotation").getElement("xsd:documentation").getText;
};

~catHelp = { arg what;
	h.getDocumentElement.getChildNodes.detect({|item| item.getAttribute("name")== (what ++ "Type")}).getElement("xsd:annotation").getLastChild.getText;
};

~catPost = {
h.getDocumentElement.getChildNodes.do({|item| item.getAttribute("name").postln;})
};
)

~catPost.value;
~catHelp.("chem_comp");
~help.("pdbx_stereo_config")
