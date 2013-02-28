// Load the xml Document
d = DOMDocument("/Users/stratosbichakis/Documents/SuperCollider files/BIOCHEMICAL-molecule-Project/03.Data/3SPX.xml");

//Create a Dictionary consisting of Lists fot atom_site data

(
var names, category;
category = d.getDocumentElement.getElement("PDBx:atom_siteCategory");
names = category.getFirstChild.getChildNodes.collect(_.getNodeName);
~atom_site = ();
~atom_site.put('id', category.getChildNodes.collectAs({|item| item.getAttribute("id").asInteger}, List));
names.do({|name| ~atom_site.put(name.asSymbol, category.getChildNodes.collectAs({|item| item.getElement(name).getText}, List))});
)

//access a specific entry, with index [i]
~atom_site.keys.do({|name| postf("%\t\t%\n", name, ~atom_site[name][4])})

// Simple SynthDef for multitude assignment
Ndef(\test, {|mul = 0.3| Out.ar(0, Pan2.ar(WhiteNoise.ar(mul), 0));})
SynthDef(\test, {|mul = 0.3| Out.ar(0, Pan2.ar(WhiteNoise.ar(mul), 0));}).add

(	//a simple Task
t = Task({
	x = Ndef(\test);
	~atom_site['id'].do({|item, i|
		x.set(\mul, ~atom_site['PDBx:B_iso_or_equiv'][i].asFloat/~atom_site['PDBx:B_iso_or_equiv'].asFloat.maxItem);
	0.1.wait;
	});
});
)

t.start; x.play

t.stop;x.free;

( 	// a simple Pattern
var list;
list = ~atom_site['PDBx:B_iso_or_equiv'].collect(_.asFloat/~atom_site['PDBx:B_iso_or_equiv'].asFloat.maxItem);
p = Pmono(\test,
	\mul, Pseq(list),
	\dur, 0.02).play
)


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
	t = ((g/abs(g))*acos(g)*180/pi).round(0.1);
};
)

(
~torsionAngle = {|seq_id|
var indexCAp, indexCp, indexN, indexCA, indexC, indexNa;		//indexes
var carCAp, carCp, carN, carCA, carC, carNa;				//Cartesian
var phi, psi, omega;

indexCAp = ~atom_site['PDBx:label_seq_id'].detectIndex({|item, i| 					~atom_site['PDBx:label_seq_id'][i].asInteger == (seq_id - 1) &&
	~atom_site['PDBx:label_atom_id'][i]=="CA"});

indexCp = ~atom_site['PDBx:label_seq_id'].detectIndex({|item, i| 					~atom_site['PDBx:label_seq_id'][i].asInteger == (seq_id - 1) &&
	~atom_site['PDBx:label_atom_id'][i]=="C"});

indexN = ~atom_site['PDBx:label_seq_id'].detectIndex({|item, i| 					~atom_site['PDBx:label_seq_id'][i].asInteger == seq_id &&
	~atom_site['PDBx:label_atom_id'][i]=="N"});

indexCA = ~atom_site['PDBx:label_seq_id'].detectIndex({|item, i| 					~atom_site['PDBx:label_seq_id'][i].asInteger == seq_id &&
	~atom_site['PDBx:label_atom_id'][i]=="CA"});

indexC = ~atom_site['PDBx:label_seq_id'].detectIndex({|item, i| 					~atom_site['PDBx:label_seq_id'][i].asInteger == seq_id &&
	~atom_site['PDBx:label_atom_id'][i]=="C"});

indexNa = ~atom_site['PDBx:label_seq_id'].detectIndex({|item, i| 					~atom_site['PDBx:label_seq_id'][i].asInteger == (seq_id + 1) &&
	~atom_site['PDBx:label_atom_id'][i]=="N"});

carCAp = Cartesian(
	x: ~atom_site['PDBx:Cartn_x'][indexCAp].asFloat,
	y: ~atom_site['PDBx:Cartn_y'][indexCAp].asFloat,
	z: ~atom_site['PDBx:Cartn_z'][indexCAp].asFloat);

carCp = Cartesian(
	x: ~atom_site['PDBx:Cartn_x'][indexCp].asFloat,
	y: ~atom_site['PDBx:Cartn_y'][indexCp].asFloat,
	z: ~atom_site['PDBx:Cartn_z'][indexCp].asFloat);

carN = Cartesian(
	x: ~atom_site['PDBx:Cartn_x'][indexN].asFloat,
	y: ~atom_site['PDBx:Cartn_y'][indexN].asFloat,
	z: ~atom_site['PDBx:Cartn_z'][indexN].asFloat);

carCA = Cartesian(
	x: ~atom_site['PDBx:Cartn_x'][indexCA].asFloat,
	y: ~atom_site['PDBx:Cartn_y'][indexCA].asFloat,
	z: ~atom_site['PDBx:Cartn_z'][indexCA].asFloat);

carC = Cartesian(
	x: ~atom_site['PDBx:Cartn_x'][indexC].asFloat,
	y: ~atom_site['PDBx:Cartn_y'][indexC].asFloat,
	z: ~atom_site['PDBx:Cartn_z'][indexC].asFloat);

carNa = Cartesian(
	x: ~atom_site['PDBx:Cartn_x'][indexNa].asFloat,
	y: ~atom_site['PDBx:Cartn_y'][indexNa].asFloat,
	z: ~atom_site['PDBx:Cartn_z'][indexNa].asFloat);

phi = ~torFromCart.value(carCp, carN, carCA, carC);
psi = ~torFromCart.value(carN, carCA, carC, carNa);
omega = ~torFromCart.value(carCAp, carCp, carN, carCA);
[phi, psi, omega];
};
)

//number of amino acids
~atom_site['PDBx:label_seq_id'].select({|item, i| item != nil}).asInteger.maxItem

318.do({|i| ~torsionAngle.(i+3).postln})
318.collect({|i| ~torsionAngle.(i+3)[0]}).plot
318.collect({|i| ~torsionAngle.(i+3)[1]}).plot
318.collect({|i| ~torsionAngle.(i+3)[2]}).plot

// isotropic displacement mean for each Amino Acid
318.collect({|seq_id| ~atom_site['PDBx:B_iso_or_equiv'].select({|item, i| ~atom_site['PDBx:label_seq_id'][i] ==(seq_id+3).asString }).collect(_.asFloat).mean.round(0.01)}).plot;