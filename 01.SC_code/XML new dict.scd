// Load the xml Document
d = DOMDocument("../03.Data/3SPX.xml".resolveRelative);


( //Create a Dictionary consisting of Ordered Lists for each atom_site data type
var names, category;
category = d.getDocumentElement.getElement("PDBx:atom_siteCategory");
names = category.getFirstChild.getChildNodes.collect(_.getNodeName);

~atom_site = ();
names.do({|name| ~atom_site.put(
		name.select({|item, i| i >= 5}).asSymbol,
		Order.newFromIndices(
			category.getChildNodes.collect({|item|
				item.getElement(name).getText}),
			category.getChildNodes.collect({|item|
				item.getAttribute("id").asInteger})))});

~atom_site.keysValuesChange({|key, elem|
	switch
	(key,
		'label_seq_id',		{elem.collect({|item| if(item!=nil,{item.asInteger})})},
		'auth_seq_id',		{elem.collect({|item| if(item!=nil,{item.asInteger})})},
		'label_entity_id',	{elem.collect({|item| if(item!=nil,{item.asInteger})})},
		'Cartn_x',			{elem.collect({|item| if(item!=nil,{item.asFloat})})},
		'Cartn_y',			{elem.collect({|item| if(item!=nil,{item.asFloat})})},
		'Cartn_z',			{elem.collect({|item| if(item!=nil,{item.asFloat})})},
		'B_iso_or_equiv',	{elem.collect({|item| if(item!=nil,{item.asFloat})})},
		'occupancy',		{elem.collect({|item| if(item!=nil,{item.asFloat})})},
		{elem})
	});
"ok";
)

( //Create a Dictionary consisting of Ordered Lists for each amino acid
var category;
category = d.getDocumentElement.getElement("PDBx:entity_poly_seqCategory");

~entity_poly_seq = Order.newFromIndices(
			category.getChildNodes.collect({|item|
				item.getAttribute("mon_id")}),
			category.getChildNodes.collect({|item|
				item.getAttribute("num").asInteger}));
"ok";
)

//access a specific entry, with index [i]
~atom_site.keys.do({|name| postf("%\t\t%\n", name, ~atom_site[name][1])});

// Calculate torsion angles phi[C-1, N, CA, C] and psi[N, CA, C, N2]
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
	t = (sign(g)*acos(g)*180/pi).round(0.1);
};
)

(
~torsionAngles = {|seq_id|
var indexCAp, indexCp, indexN, indexCA, indexC, indexNa;		//indeces
var carCAp, carCp, carN, carCA, carC, carNa;				//Cartesian
var phi, psi, omega;

indexCAp = ~atom_site['label_seq_id'].detectIndex({|item, i|
		item == (seq_id - 1) &&
		~atom_site['label_atom_id'][i]=="CA"});

indexCp = ~atom_site['label_seq_id'].detectIndex({|item, i|
		item == (seq_id - 1) &&
		~atom_site['label_atom_id'][i]=="C"});

indexN = ~atom_site['label_seq_id'].detectIndex({|item, i|
		item == seq_id &&
		~atom_site['label_atom_id'][i]=="N"});

indexCA = ~atom_site['label_seq_id'].detectIndex({|item, i|
		item == seq_id &&
		~atom_site['label_atom_id'][i]=="CA"});

indexC = ~atom_site['label_seq_id'].detectIndex({|item, i|
		item == seq_id &&
		~atom_site['label_atom_id'][i]=="C"});

indexNa = ~atom_site['label_seq_id'].detectIndex({|item, i|
		item == (seq_id + 1) &&
		~atom_site['label_atom_id'][i]=="N"});

carCAp = Cartesian(
	x: ~atom_site['Cartn_x'][indexCAp],
	y: ~atom_site['Cartn_y'][indexCAp],
	z: ~atom_site['Cartn_z'][indexCAp]);

carCp = Cartesian(
	x: ~atom_site['Cartn_x'][indexCp],
	y: ~atom_site['Cartn_y'][indexCp],
	z: ~atom_site['Cartn_z'][indexCp]);

carN = Cartesian(
	x: ~atom_site['Cartn_x'][indexN],
	y: ~atom_site['Cartn_y'][indexN],
	z: ~atom_site['Cartn_z'][indexN]);

carCA = Cartesian(
	x: ~atom_site['Cartn_x'][indexCA],
	y: ~atom_site['Cartn_y'][indexCA],
	z: ~atom_site['Cartn_z'][indexCA]);

carC = Cartesian(
	x: ~atom_site['Cartn_x'][indexC],
	y: ~atom_site['Cartn_y'][indexC],
	z: ~atom_site['Cartn_z'][indexC]);

carNa = Cartesian(
	x: ~atom_site['Cartn_x'][indexNa],
	y: ~atom_site['Cartn_y'][indexNa],
	z: ~atom_site['Cartn_z'][indexNa]);

phi = ~torFromCart.value(carCp, carN, carCA, carC);
psi = ~torFromCart.value(carN, carCA, carC, carNa);
omega = ~torFromCart.value(carCAp, carCp, carN, carCA);
[phi, psi, omega];
};
)

//number of amino acids
~atom_site['label_seq_id'].select({|item| item != nil}).maxItem
//or
~atom_site['label_seq_id'].at(~atom_site['label_seq_id'].lastIndexForWhich({|item| item != nil}))

~torsionAnglesList = Order.fill(318, {|i| ~torsionAngles.(i+3)})
~torsionAnglesList.flop.plot


// isotropic displacement mean for each Amino Acid
~b_isoList = 318.collect({|seq_id| ~atom_site['B_iso_or_equiv'].select({|item, i| ~atom_site['label_seq_id'][i] ==(seq_id+3)}).mean.round(0.01)});



(	// Included ligands
var category;
category = d.getDocumentElement.getElement("PDBx:chem_compCategory");
~chem_comps = category.getChildNodes.collect({|item|
				item.getAttribute("id")});
"ok";
)

// Access the specific xml files for each amino Acid
(
~chem_compsDocuments = Dictionary.fill(~chem_comps.size,
	{|i| ~chem_comps[i].asSymbol ->
		DOMDocument("~/Downloads/components-pub-xml/".standardizePath +/+
			~chem_comps[i] ++".xml")} );
)
~chem_compsDocuments[\ALA]


~sheetCat = d.getElementsByTagName("PDBx:struct_sheet_range");
~sheetBeg = ~sheetCat.collect({|item, i| item.getElement 						("PDBx:beg_label_seq_id").getText.asInteger;});
~sheetEnd = ~sheetCat.collect({|item, i| item.getElement 						("PDBx:end_label_seq_id").getText.asInteger;});
~helixCat = d.getElementsByTagName("PDBx:struct_conf");
~helixBeg = ~helixCat.collect({|item, i| item.getElement 						("PDBx:beg_label_seq_id").getText.asInteger;});
~helixEnd = ~helixCat.collect({|item, i| item.getElement 						("PDBx:end_label_seq_id").getText.asInteger;});



//-------------------------Sound Examples---------------------------
// Simple SynthDef for multitude assignment
Ndef(\test, {|mul = 0.3| Out.ar(0, Pan2.ar(WhiteNoise.ar(mul), 0));})
SynthDef(\test, {|mul = 0.3| Out.ar(0, Pan2.ar(WhiteNoise.ar(mul), 0));}).add
(
SynthDef(\sine, {|mul = 0.58, angle = 0, fb, pan = 0|
	Out.ar(0, Pan2.ar(HPF.ar(SinOsc.ar(45+(angle*0.2), 0, mul: mul)*SinOsc.ar([85, 110]).clip(-0.6, 0.6), 1000), pan));}, [0, 1, 1]).add
)
(//a simple Task
x = Ndef(\test);
t = Task({loop{
	~atom_site['B_iso_or_equiv'].doRange({|item, i|
		x.set(\mul, item/~atom_site['B_iso_or_equiv'].maxItem);
	0.1.wait;
		},
		from: ~atom_site['label_seq_id'].detectIndex({|item| item == 2}),
	      to: ~atom_site['label_seq_id'].lastIndexForWhich({|item| item == 22}));
}});
)

(//a simple Task of angles Sonification

t = Task({
	x = Array.fill(3, {Synth(\sine)});

	~torsionAnglesList.doRange({|item, i|
			x[0].set(\angle, item[0], \fb, ~b_isoList[i], \pan, -1);
			x[1].set(\angle, item[1], \fb, ~b_isoList[i], \pan, 1);
			x[2].set(\angle, item[2], \fb, ~b_isoList[i], \pan, 0, \mul, 0);
			// To Do: add percussive sounds on helix and sheet beginning

	0.12.wait;
}, 0,310 );
1.wait;
	x.do(_.free);
});
)

t.start;
t.reset;
x.do(_.free);

{Silent.ar()!.scope

( 	// a simple Pattern
var list;
list = ~atom_site['B_iso_or_equiv'].collectAs(_/~atom_site['B_iso_or_equiv'].maxItem, List);
p = Pmono(\test,
	\mul, Pseq(list),
	\dur, 0.02).play
)
