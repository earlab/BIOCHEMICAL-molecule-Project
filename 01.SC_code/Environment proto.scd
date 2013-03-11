q = Environment.make({|env|
	~init = {
		env.construct_atom_site;
		env.construct_torsionAnglesList;
		env.construct_b_isoList;
	};

	~pdbml = DOMDocument("../03.Data/3SPX.xml".resolveRelative);

	~construct_atom_site = {|env|
		var names, category;
		category = env.pdbml.getDocumentElement.getElement("PDBx:atom_siteCategory");
		names = category.getFirstChild.getChildNodes.collect(_.getNodeName);

		env.atom_site = ();
		names.do({|name| env.atom_site.put(
			name.drop(5).asSymbol,
				Order.newFromIndices(
					category.getChildNodes.collect({|item|
						item.getElement(name).getText}),
					category.getChildNodes.collect({|item|
						item.getAttribute("id").asInteger})))});

env.atom_site.keysValuesChange({|key, elem|
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
	};

	~torsionAngles = {|env, seq_id|
var indexCAp, indexCp, indexN, indexCA, indexC, indexNa;		//indeces
var carCAp, carCp, carN, carCA, carC, carNa;				//Cartesian
var phi, psi, omega;

indexCAp = env.atom_site['label_seq_id'].detectIndex({|item, i|
		item == (seq_id - 1) &&
		env.atom_site['label_atom_id'][i]=="CA"});

indexCp = env.atom_site['label_seq_id'].detectIndex({|item, i|
		item == (seq_id - 1) &&
		env.atom_site['label_atom_id'][i]=="C"});

indexN = env.atom_site['label_seq_id'].detectIndex({|item, i|
		item == seq_id &&
		env.atom_site['label_atom_id'][i]=="N"});

indexCA = env.atom_site['label_seq_id'].detectIndex({|item, i|
		item == seq_id &&
		env.atom_site['label_atom_id'][i]=="CA"});

indexC = env.atom_site['label_seq_id'].detectIndex({|item, i|
		item == seq_id &&
		env.atom_site['label_atom_id'][i]=="C"});

indexNa = env.atom_site['label_seq_id'].detectIndex({|item, i|
		item == (seq_id + 1) &&
		env.atom_site['label_atom_id'][i]=="N"});

carCAp = Cartesian(
	x: env.atom_site['Cartn_x'][indexCAp],
	y: env.atom_site['Cartn_y'][indexCAp],
	z: env.atom_site['Cartn_z'][indexCAp]);

carCp = Cartesian(
	x: env.atom_site['Cartn_x'][indexCp],
	y: env.atom_site['Cartn_y'][indexCp],
	z: env.atom_site['Cartn_z'][indexCp]);

carN = Cartesian(
	x: env.atom_site['Cartn_x'][indexN],
	y: env.atom_site['Cartn_y'][indexN],
	z: env.atom_site['Cartn_z'][indexN]);

carCA = Cartesian(
	x: env.atom_site['Cartn_x'][indexCA],
	y: env.atom_site['Cartn_y'][indexCA],
	z: env.atom_site['Cartn_z'][indexCA]);

carC = Cartesian(
	x: env.atom_site['Cartn_x'][indexC],
	y: env.atom_site['Cartn_y'][indexC],
	z: env.atom_site['Cartn_z'][indexC]);

carNa = Cartesian(
	x: env.atom_site['Cartn_x'][indexNa],
	y: env.atom_site['Cartn_y'][indexNa],
	z: env.atom_site['Cartn_z'][indexNa]);

phi = env.torFromCart(carCp, carN, carCA, carC);
psi = env.torFromCart(carN, carCA, carC, carNa);
omega = env.torFromCart(carCAp, carCp, carN, carCA);
[phi, psi, omega];
};

	~angleFromCart = {|env, a, b, c|
		acos(((a-b)*(c-b)).asArray.sum/(dist(a,b)*dist(b,c)))*180/pi};

	~crossFromCart = {|env, a, b, c|
		var e, f;
		e = (b-a).asArray *.t (c-b).asArray;
		f = [e[1][2]-e[2][1], e[2][0]-e[0][2], e[0][1]-e[1][0]];};

	~torFromCart = {|env, a, b, c, d|
		var e,f,g,t;
		e = env.crossFromCart(a, b, c);
		f = env.crossFromCart(d, c, b);
		g = (e*f).sum/(sqrt(e.collect(_.squared).sum)*sqrt(f.collect(_.squared).sum));
		t = (sign(g)*acos(g)*180/pi).round(0.1);};

	~construct_torsionAnglesList = {env.torsionAnglesList = Order.fill(318, {|i| env.torsionAngles(i+3)});};

	~construct_b_isoList = {env.b_isoList = 318.collect({|seq_id| env.atom_site['B_iso_or_equiv'].select({|item, i| env.atom_site['label_seq_id'][i] ==(seq_id+3)}).mean.round(0.01)});};


}).know_(true);

q.init

Archive.global.put(\myData, q);