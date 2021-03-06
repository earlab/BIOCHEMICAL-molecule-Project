#include "AtomController.h"

AtomController::AtomController() {

}

void AtomController::addAtom(Atom atom) {
	atoms.push_back(atom);
}

void AtomController::draw() {
	lastAtomPosition = ofVec3f(0,0,0);  //we use these to draw lines
	lastAtomGroup = 0;					//between atoms belonging to the same aminoacid
	for (list<Atom>::iterator atom = atoms.begin(); atom != atoms.end(); atom++){
		if (atom->tempTransparency > atom->transparency){		//Newly created atoms start out more opaque
			atom->tempTransparency = atom->tempTransparency - 5;//and fade out till they reach their normal transparency
		}
		ofSetColor(atom->color.r, atom->color.g,atom->color.b, atom->tempTransparency);
		ofVec3f position = atom->position;

		position.x = position.x + cos(ofGetElapsedTimef()*sqrt(2*atom->displacement));
		position.y = position.y + sin(ofGetElapsedTimef()*sqrt(2*atom->displacement));
		position.z = position.z + 0.5*cos(-ofGetElapsedTimef()*sqrt(2*atom->displacement)) + 0.5*sin(-ofGetElapsedTimef()*sqrt(2*atom->displacement));

		ofSphere(position, atom->displacement*2); 
		if (atom != atoms.begin() && atom->group == lastAtomGroup) {
			ofSetColor(128,255,0,96);  //Connect atoms belonging to the same
			ofLine(lastAtomPosition,atom->position);
		}

		if (atom != atoms.begin() && atom->group != lastAtomGroup) {
			ofSetColor(231,245,152,96);
			ofSetLineWidth(10.0);
			ofLine(lastAcidPosition, atom->position);
			ofSetColor(231,245,152,60);
			ofSphere(atom->position,1.5);
			ofSetLineWidth(1.0);
		}

		if (atom->group != lastAtomGroup) {
			lastAcidPosition = atom->position;
		}

		lastAtomPosition = atom->position;
		lastAtomGroup = atom->group;
	}
}

void AtomController::clear() { //we erase all atoms created so far
	for (list<Atom>::iterator atom = atoms.begin(); atom != atoms.end(); atom++) {
		atom = atoms.erase(atom);
	}
}