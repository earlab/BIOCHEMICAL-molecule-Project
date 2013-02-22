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

//access a specific entry
~atom_site.keys.do({|name| postf("%\t\t%\n", name, ~atom_site[name][0])})

// Simple SynthDef for multitude assignment
Ndef(\test, {|mul = 0.3| Out.ar(0, Pan2.ar(WhiteNoise.ar(mul), 0));})
SynthDef(\test, {|mul = 0.3| Out.ar(0, Pan2.ar(WhiteNoise.ar(mul), 0));}).add

x = Ndef(\test)

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